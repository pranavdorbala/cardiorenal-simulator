"""
Checkpoint management for HPC training with Slurm preemption support.

Handles:
  - Automatic checkpoint saving/loading with versioning
  - Slurm SIGTERM/SIGUSR1 signal handling for graceful preemption
  - Training state preservation (model, optimizer, scheduler, epoch, RNG states)
  - Checkpoint cleanup (keep N most recent + best)
"""

import os
import signal
import sys
import glob
import json
import logging
import torch
import numpy as np
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)


class CheckpointManager:
    """Manages training checkpoints with Slurm-aware signal handling."""

    def __init__(self, save_dir, keep_last_n=3, keep_best=True):
        """
        Args:
            save_dir: Directory for checkpoint files.
            keep_last_n: Number of recent periodic checkpoints to retain.
            keep_best: Whether to keep the best-validation checkpoint.
        """
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)
        self.keep_last_n = keep_last_n
        self.keep_best = keep_best

        self._preempt_flag = False
        self._original_sigterm = None
        self._original_sigusr1 = None

    # ------------------------------------------------------------------
    # Signal handling for Slurm preemption
    # ------------------------------------------------------------------

    def register_signal_handlers(self):
        """Install handlers for SIGTERM and SIGUSR1 (Slurm preemption signals).

        Slurm sends SIGUSR1 (configurable) or SIGTERM before killing a job.
        The handler sets a flag so the training loop can save a checkpoint
        and exit cleanly.
        """
        self._original_sigterm = signal.getsignal(signal.SIGTERM)
        self._original_sigusr1 = signal.getsignal(signal.SIGUSR1)

        signal.signal(signal.SIGTERM, self._handle_preempt)
        signal.signal(signal.SIGUSR1, self._handle_preempt)
        logger.info("Registered Slurm preemption signal handlers (SIGTERM, SIGUSR1)")

    def _handle_preempt(self, signum, frame):
        """Set the preemption flag; training loop checks this each epoch."""
        sig_name = signal.Signals(signum).name
        logger.warning(f"Received {sig_name} — flagging for checkpoint and exit")
        self._preempt_flag = True

    @property
    def should_preempt(self):
        """Check if a preemption signal has been received."""
        return self._preempt_flag

    def restore_signal_handlers(self):
        """Restore original signal handlers."""
        if self._original_sigterm is not None:
            signal.signal(signal.SIGTERM, self._original_sigterm)
        if self._original_sigusr1 is not None:
            signal.signal(signal.SIGUSR1, self._original_sigusr1)

    # ------------------------------------------------------------------
    # Save / Load
    # ------------------------------------------------------------------

    def save_checkpoint(self, state, epoch, val_loss, is_best=False, tag=None):
        """Save a training checkpoint.

        Args:
            state: Dict with model_state_dict, optimizer_state_dict,
                   scheduler_state_dict, scaler_state_dict, etc.
            epoch: Current epoch number.
            val_loss: Validation loss for this epoch.
            is_best: If True, also save as best_model.pt.
            tag: Optional tag (e.g., 'preempt', 'final').
        """
        # Build filename
        if tag:
            filename = f"checkpoint_{tag}.pt"
        else:
            filename = f"checkpoint_epoch{epoch:04d}.pt"
        filepath = self.save_dir / filename

        # Augment state with metadata
        state['epoch'] = epoch
        state['val_loss'] = val_loss
        state['timestamp'] = datetime.now().isoformat()

        # Save RNG states for exact reproducibility on resume
        state['rng_python'] = torch.random.get_rng_state()
        state['rng_numpy'] = np.random.get_state()
        if torch.cuda.is_available():
            state['rng_cuda'] = torch.cuda.get_rng_state_all()

        torch.save(state, filepath)
        logger.info(f"Saved checkpoint: {filepath} (epoch={epoch}, val_loss={val_loss:.6f})")

        if is_best and self.keep_best:
            best_path = self.save_dir / "best_model.pt"
            torch.save(state, best_path)
            logger.info(f"Saved best model: {best_path}")

        # Cleanup old periodic checkpoints (not best, not tagged)
        self._cleanup_old_checkpoints()

        return filepath

    def _cleanup_old_checkpoints(self):
        """Remove old periodic checkpoints, keeping only the last N."""
        pattern = str(self.save_dir / "checkpoint_epoch*.pt")
        ckpts = sorted(glob.glob(pattern))
        if len(ckpts) > self.keep_last_n:
            for old in ckpts[:-self.keep_last_n]:
                os.remove(old)
                logger.debug(f"Removed old checkpoint: {old}")

    def load_checkpoint(self, path=None, map_location=None):
        """Load a checkpoint. If path is None, loads the latest available.

        Priority: preempt checkpoint > latest epoch checkpoint > best model.

        Returns:
            Checkpoint dict, or None if no checkpoint found.
        """
        if path is not None:
            if os.path.exists(path):
                logger.info(f"Loading specified checkpoint: {path}")
                return torch.load(path, map_location=map_location, weights_only=False)
            else:
                logger.warning(f"Specified checkpoint not found: {path}")
                return None

        # Auto-detect: preempt > latest epoch > best
        candidates = []

        preempt_path = self.save_dir / "checkpoint_preempt.pt"
        if preempt_path.exists():
            candidates.append(("preempt", preempt_path))

        epoch_pattern = str(self.save_dir / "checkpoint_epoch*.pt")
        epoch_ckpts = sorted(glob.glob(epoch_pattern))
        if epoch_ckpts:
            candidates.append(("latest_epoch", Path(epoch_ckpts[-1])))

        best_path = self.save_dir / "best_model.pt"
        if best_path.exists():
            candidates.append(("best", best_path))

        if not candidates:
            logger.info("No checkpoint found — starting from scratch")
            return None

        label, chosen = candidates[0]
        logger.info(f"Resuming from {label} checkpoint: {chosen}")
        return torch.load(chosen, map_location=map_location, weights_only=False)

    def get_resume_epoch(self, checkpoint):
        """Extract the epoch to resume from (next epoch after saved)."""
        if checkpoint is None:
            return 1
        return checkpoint.get('epoch', 0) + 1

    # ------------------------------------------------------------------
    # Restore RNG states
    # ------------------------------------------------------------------

    def restore_rng_states(self, checkpoint):
        """Restore Python/NumPy/CUDA RNG states from checkpoint."""
        if checkpoint is None:
            return
        if 'rng_python' in checkpoint:
            torch.random.set_rng_state(checkpoint['rng_python'])
        if 'rng_numpy' in checkpoint:
            np.random.set_state(checkpoint['rng_numpy'])
        if 'rng_cuda' in checkpoint and torch.cuda.is_available():
            torch.cuda.set_rng_state_all(checkpoint['rng_cuda'])
        logger.info("Restored RNG states from checkpoint")

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def list_checkpoints(self):
        """List all checkpoints in save_dir."""
        all_ckpts = sorted(glob.glob(str(self.save_dir / "*.pt")))
        info = []
        for p in all_ckpts:
            try:
                ckpt = torch.load(p, map_location='cpu', weights_only=False)
                info.append({
                    'path': p,
                    'epoch': ckpt.get('epoch', '?'),
                    'val_loss': ckpt.get('val_loss', '?'),
                    'timestamp': ckpt.get('timestamp', '?'),
                })
            except Exception:
                info.append({'path': p, 'error': 'could not load'})
        return info

    def save_preempt_checkpoint(self, state, epoch, val_loss):
        """Save a checkpoint specifically tagged for preemption resume."""
        return self.save_checkpoint(state, epoch, val_loss, tag='preempt')
