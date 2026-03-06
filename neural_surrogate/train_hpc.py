"""
HPC training script for the Neural ODE cardiorenal surrogate.

Supports:
  - PyTorch Distributed Data Parallel (DDP) for multi-GPU / multi-node
  - Mixed precision training (float32 ODE solve, float16 where safe)
  - Slurm preemption handling with automatic checkpoint/resume
  - Curriculum learning (short windows -> long windows)
  - Gradient accumulation for effective larger batch sizes
  - TensorBoard logging
  - NaN/Inf detection with automatic recovery

Usage:
    # Single GPU
    python -m neural_surrogate.train_hpc --data training_data.h5

    # Multi-GPU (launched by Slurm or torchrun)
    torchrun --nproc_per_node=4 -m neural_surrogate.train_hpc --data training_data.h5

    # Resume from preempted job
    python -m neural_surrogate.train_hpc --data training_data.h5 --resume auto
"""

import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import DataLoader
from torch.utils.data.distributed import DistributedSampler

import numpy as np
import h5py
import os
import sys
import time
import argparse
import json
import logging
import math
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from neural_surrogate.model import (
    HallowSurrogate, CLINICAL_KEYS, CLINICAL_WEIGHTS,
    N_CLINICAL, N_CRITICAL, CRITICAL_STATE_INDICES,
    LATENT_DIM, N_STATE, N_KNOBS,
)
from neural_surrogate.train import HallowTrajectoryDataset, build_weight_tensor
from neural_surrogate.checkpoint_manager import CheckpointManager

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(rank, log_dir):
    """Configure logging — only rank 0 logs to console + file."""
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger()
    logger.setLevel(logging.INFO if rank == 0 else logging.WARNING)

    fmt = logging.Formatter(
        f"[%(asctime)s][Rank {rank}][%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if rank == 0:
        # Console
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(fmt)
        logger.addHandler(ch)

        # File
        fh = logging.FileHandler(log_dir / "train.log", mode='a')
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


# ---------------------------------------------------------------------------
# DDP setup
# ---------------------------------------------------------------------------

def setup_distributed():
    """Initialize DDP if running under torchrun / Slurm.

    Returns:
        (rank, local_rank, world_size, device)
    """
    if 'RANK' in os.environ:
        # Launched via torchrun or Slurm srun with torchrun
        rank = int(os.environ['RANK'])
        local_rank = int(os.environ['LOCAL_RANK'])
        world_size = int(os.environ['WORLD_SIZE'])
    elif 'SLURM_PROCID' in os.environ:
        # Launched via srun directly (no torchrun)
        rank = int(os.environ['SLURM_PROCID'])
        local_rank = int(os.environ.get('SLURM_LOCALID', rank))
        world_size = int(os.environ['SLURM_NTASKS'])
        # Set environment for torch.distributed
        os.environ['RANK'] = str(rank)
        os.environ['LOCAL_RANK'] = str(local_rank)
        os.environ['WORLD_SIZE'] = str(world_size)
        if 'MASTER_ADDR' not in os.environ:
            import subprocess
            nodelist = os.environ.get('SLURM_NODELIST', 'localhost')
            # Resolve first node
            try:
                result = subprocess.run(
                    ['scontrol', 'show', 'hostnames', nodelist],
                    capture_output=True, text=True
                )
                master = result.stdout.strip().split('\n')[0]
            except Exception:
                master = 'localhost'
            os.environ['MASTER_ADDR'] = master
        if 'MASTER_PORT' not in os.environ:
            os.environ['MASTER_PORT'] = str(29500 + int(os.environ.get('SLURM_JOBID', '0')) % 1000)
    else:
        # Single-process
        rank = 0
        local_rank = 0
        world_size = 1

    if world_size > 1:
        backend = 'nccl' if torch.cuda.is_available() else 'gloo'
        dist.init_process_group(backend=backend, rank=rank, world_size=world_size)
        if torch.cuda.is_available():
            torch.cuda.set_device(local_rank)
            device = torch.device(f'cuda:{local_rank}')
        else:
            device = torch.device('cpu')
    else:
        if torch.cuda.is_available():
            device = torch.device('cuda:0')
        else:
            device = torch.device('cpu')

    return rank, local_rank, world_size, device


def cleanup_distributed():
    if dist.is_initialized():
        dist.destroy_process_group()


# ---------------------------------------------------------------------------
# Curriculum schedule
# ---------------------------------------------------------------------------

def get_curriculum_window(epoch, total_epochs, min_window=8, max_window=48):
    """Linearly increase window size over training.

    Start with short windows (easier for ODE solver, faster convergence)
    and gradually lengthen to capture long-range dynamics.
    """
    warmup_frac = 0.3  # spend first 30% of training ramping up
    progress = min(epoch / (total_epochs * warmup_frac), 1.0)
    window = int(min_window + progress * (max_window - min_window))
    # Round to nearest multiple of 4 for clean batching
    window = max(min_window, (window // 4) * 4)
    return min(window, max_window)


# ---------------------------------------------------------------------------
# Training step
# ---------------------------------------------------------------------------

def train_one_epoch(model, loader, optimizer, clinical_weights, device,
                    grad_accum_steps=1, max_grad_norm=10.0):
    """Run one training epoch. Returns average loss."""
    model.train()
    total_loss = 0.0
    n_batches = 0
    nan_batches = 0

    optimizer.zero_grad()

    for step, batch in enumerate(loader):
        state = batch['state'].to(device, non_blocking=True)
        knobs = batch['knobs'].to(device, non_blocking=True)
        t_eval = batch['times'][0].to(device, non_blocking=True)
        clinical_gt = batch['clinical'].to(device, non_blocking=True)
        crit_state_gt = batch['critical_state'].to(device, non_blocking=True)

        # Forward
        try:
            clinical_pred, crit_pred = model(state, knobs, t_eval)
        except RuntimeError as e:
            if 'nan' in str(e).lower() or 'inf' in str(e).lower():
                nan_batches += 1
                logging.warning(f"NaN/Inf in forward pass (batch {step}), skipping")
                continue
            raise

        # odeint returns (T, B, ...), permute to (B, T, ...)
        clinical_pred = clinical_pred.permute(1, 0, 2)
        crit_pred = crit_pred.permute(1, 0, 2)

        # NaN check on predictions
        if torch.isnan(clinical_pred).any() or torch.isnan(crit_pred).any():
            nan_batches += 1
            logging.warning(f"NaN in predictions (batch {step}), skipping")
            optimizer.zero_grad()
            continue

        # Weighted clinical loss
        clinical_err = (clinical_pred - clinical_gt) ** 2
        clinical_loss = (clinical_err * clinical_weights.unsqueeze(0).unsqueeze(0)).mean()

        # Critical state loss
        crit_loss = ((crit_pred - crit_state_gt) ** 2).mean()

        # ODE regularization
        unwrapped = model.module if hasattr(model, 'module') else model
        z0 = unwrapped.encoder(state, knobs)
        unwrapped.ode_func.set_knobs(knobs)
        dzdt = unwrapped.ode_func(torch.tensor(0.0, device=device), z0)
        reg_loss = 0.01 * (dzdt ** 2).mean()

        loss = clinical_loss + 0.5 * crit_loss + reg_loss

        # Scale for gradient accumulation
        loss_scaled = loss / grad_accum_steps
        loss_scaled.backward()

        if (step + 1) % grad_accum_steps == 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_grad_norm)
            optimizer.step()
            optimizer.zero_grad()

        total_loss += loss.item()
        n_batches += 1

    # Handle remaining gradients
    if (step + 1) % grad_accum_steps != 0:
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_grad_norm)
        optimizer.step()
        optimizer.zero_grad()

    avg_loss = total_loss / max(n_batches, 1)

    if nan_batches > 0:
        logging.warning(f"Skipped {nan_batches} NaN batches this epoch")

    return avg_loss, nan_batches


@torch.no_grad()
def validate(model, loader, clinical_weights, device):
    """Run validation. Returns average loss and per-output losses."""
    model.eval()
    total_loss = 0.0
    per_output_loss = torch.zeros(N_CLINICAL, device=device)
    n_batches = 0

    for batch in loader:
        state = batch['state'].to(device, non_blocking=True)
        knobs = batch['knobs'].to(device, non_blocking=True)
        t_eval = batch['times'][0].to(device, non_blocking=True)
        clinical_gt = batch['clinical'].to(device, non_blocking=True)
        crit_state_gt = batch['critical_state'].to(device, non_blocking=True)

        try:
            clinical_pred, crit_pred = model(state, knobs, t_eval)
        except RuntimeError:
            continue

        clinical_pred = clinical_pred.permute(1, 0, 2)
        crit_pred = crit_pred.permute(1, 0, 2)

        if torch.isnan(clinical_pred).any():
            continue

        clinical_err = (clinical_pred - clinical_gt) ** 2
        clinical_loss = (clinical_err * clinical_weights.unsqueeze(0).unsqueeze(0)).mean()
        crit_loss = ((crit_pred - crit_state_gt) ** 2).mean()
        loss = clinical_loss + 0.5 * crit_loss

        # Track per-output MSE (unweighted)
        per_output_loss += clinical_err.mean(dim=(0, 1))

        total_loss += loss.item()
        n_batches += 1

    avg_loss = total_loss / max(n_batches, 1)
    per_output_loss = per_output_loss / max(n_batches, 1)
    return avg_loss, per_output_loss


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def train_hpc(args):
    rank, local_rank, world_size, device = setup_distributed()
    is_main = (rank == 0)

    # Directories
    run_dir = Path(args.run_dir)
    ckpt_dir = run_dir / "checkpoints"
    log_dir = run_dir / "logs"

    logger = setup_logging(rank, log_dir)

    if is_main:
        logger.info("=" * 60)
        logger.info("Neural ODE Cardiorenal Surrogate — HPC Training")
        logger.info("=" * 60)
        logger.info(f"Device: {device} | World size: {world_size}")
        logger.info(f"Run directory: {run_dir}")
        logger.info(f"Args: {vars(args)}")

    # Checkpoint manager (only rank 0 saves, but all ranks need the manager for signals)
    ckpt_mgr = CheckpointManager(
        save_dir=ckpt_dir,
        keep_last_n=args.keep_checkpoints,
    )
    ckpt_mgr.register_signal_handlers()

    # TensorBoard (optional, rank 0 only)
    tb_writer = None
    if is_main:
        try:
            from torch.utils.tensorboard import SummaryWriter
            tb_writer = SummaryWriter(log_dir=str(log_dir / "tensorboard"))
            logger.info("TensorBoard logging enabled")
        except ImportError:
            logger.info("TensorBoard not available — logging to file only")

    # ---------------------------------------------------------------
    # Data
    # ---------------------------------------------------------------
    initial_window = args.window_min if args.curriculum else args.window_size
    stride = max(1, initial_window // 4)

    train_ds = HallowTrajectoryDataset(
        args.data, window_size=initial_window, stride=stride, split='train',
        train_frac=args.train_frac,
    )
    val_ds = HallowTrajectoryDataset(
        args.data, window_size=initial_window, stride=stride, split='val',
        train_frac=args.train_frac,
    )
    # Ensure val uses train normalization stats
    for attr in ('clinical_mean', 'clinical_std', 'state_mean', 'state_std',
                 'knob_mean', 'knob_std'):
        setattr(val_ds, attr, getattr(train_ds, attr))

    # Normalization stats to save with model
    norm_stats = {
        'clinical_mean': train_ds.clinical_mean.tolist(),
        'clinical_std': train_ds.clinical_std.tolist(),
        'state_mean': train_ds.state_mean.tolist(),
        'state_std': train_ds.state_std.tolist(),
        'knob_mean': train_ds.knob_mean.tolist(),
        'knob_std': train_ds.knob_std.tolist(),
    }

    # Distributed sampler
    train_sampler = DistributedSampler(train_ds, num_replicas=world_size,
                                        rank=rank, shuffle=True) if world_size > 1 else None
    val_sampler = DistributedSampler(val_ds, num_replicas=world_size,
                                      rank=rank, shuffle=False) if world_size > 1 else None

    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size,
        sampler=train_sampler,
        shuffle=(train_sampler is None),
        num_workers=args.num_workers,
        pin_memory=torch.cuda.is_available(),
        drop_last=True,
    )
    val_loader = DataLoader(
        val_ds, batch_size=args.batch_size,
        sampler=val_sampler,
        shuffle=False,
        num_workers=args.num_workers,
        pin_memory=torch.cuda.is_available(),
    )

    if is_main:
        logger.info(f"Train: {len(train_ds)} windows | Val: {len(val_ds)} windows")

    # ---------------------------------------------------------------
    # Model
    # ---------------------------------------------------------------
    model = HallowSurrogate(latent_dim=args.latent_dim).to(device)

    if world_size > 1:
        model = DDP(model, device_ids=[local_rank] if torch.cuda.is_available() else None,
                    find_unused_parameters=True)

    if is_main:
        n_params = sum(p.numel() for p in model.parameters())
        logger.info(f"Model parameters: {n_params:,}")

    # ---------------------------------------------------------------
    # Optimizer & Scheduler
    # ---------------------------------------------------------------
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=args.lr, weight_decay=args.weight_decay,
    )

    # Warmup + cosine decay
    def lr_lambda(epoch):
        if epoch < args.warmup_epochs:
            return (epoch + 1) / args.warmup_epochs
        progress = (epoch - args.warmup_epochs) / max(1, args.epochs - args.warmup_epochs)
        return 0.5 * (1 + math.cos(math.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    clinical_weights = build_weight_tensor().to(device)

    # ---------------------------------------------------------------
    # Resume from checkpoint
    # ---------------------------------------------------------------
    start_epoch = 1
    best_val_loss = float('inf')
    history = {'train_loss': [], 'val_loss': [], 'lr': [], 'window_size': []}

    if args.resume:
        resume_path = None if args.resume == 'auto' else args.resume
        ckpt = ckpt_mgr.load_checkpoint(path=resume_path, map_location=device)
        if ckpt is not None:
            unwrapped = model.module if hasattr(model, 'module') else model
            unwrapped.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            if 'scheduler_state_dict' in ckpt:
                scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt_mgr.get_resume_epoch(ckpt)
            best_val_loss = ckpt.get('best_val_loss', ckpt.get('val_loss', float('inf')))
            if 'history' in ckpt:
                history = ckpt['history']
            ckpt_mgr.restore_rng_states(ckpt)
            if is_main:
                logger.info(f"Resumed from epoch {start_epoch - 1}, "
                            f"best_val_loss={best_val_loss:.6f}")

    # ---------------------------------------------------------------
    # Training loop
    # ---------------------------------------------------------------
    if is_main:
        logger.info(f"Starting training: epochs {start_epoch}..{args.epochs}")

    consecutive_nan_epochs = 0

    for epoch in range(start_epoch, args.epochs + 1):
        epoch_t0 = time.time()

        # Curriculum: update window size and rebuild loaders if needed
        if args.curriculum:
            new_window = get_curriculum_window(
                epoch, args.epochs, args.window_min, args.window_max)
            if new_window != train_ds.window_size:
                stride = max(1, new_window // 4)
                train_ds = HallowTrajectoryDataset(
                    args.data, window_size=new_window, stride=stride,
                    split='train', train_frac=args.train_frac)
                val_ds = HallowTrajectoryDataset(
                    args.data, window_size=new_window, stride=stride,
                    split='val', train_frac=args.train_frac)
                for attr in ('clinical_mean', 'clinical_std', 'state_mean',
                             'state_std', 'knob_mean', 'knob_std'):
                    setattr(val_ds, attr, getattr(train_ds, attr))

                train_sampler = DistributedSampler(
                    train_ds, num_replicas=world_size, rank=rank,
                    shuffle=True) if world_size > 1 else None
                val_sampler = DistributedSampler(
                    val_ds, num_replicas=world_size, rank=rank,
                    shuffle=False) if world_size > 1 else None

                train_loader = DataLoader(
                    train_ds, batch_size=args.batch_size,
                    sampler=train_sampler,
                    shuffle=(train_sampler is None),
                    num_workers=args.num_workers,
                    pin_memory=torch.cuda.is_available(),
                    drop_last=True)
                val_loader = DataLoader(
                    val_ds, batch_size=args.batch_size,
                    sampler=val_sampler, shuffle=False,
                    num_workers=args.num_workers,
                    pin_memory=torch.cuda.is_available())

                if is_main:
                    logger.info(f"Curriculum: window_size={new_window}, "
                                f"train={len(train_ds)} windows")

        # Set epoch for distributed sampler
        if train_sampler is not None:
            train_sampler.set_epoch(epoch)

        # Train
        train_loss, nan_count = train_one_epoch(
            model, train_loader, optimizer, clinical_weights, device,
            grad_accum_steps=args.grad_accum, max_grad_norm=args.max_grad_norm,
        )

        # NaN recovery
        if math.isnan(train_loss) or math.isinf(train_loss):
            consecutive_nan_epochs += 1
            if is_main:
                logger.warning(f"Epoch {epoch}: NaN/Inf loss "
                               f"({consecutive_nan_epochs} consecutive)")
            if consecutive_nan_epochs >= 5:
                if is_main:
                    logger.error("5 consecutive NaN epochs — halving LR and "
                                 "reloading best checkpoint")
                ckpt = ckpt_mgr.load_checkpoint(
                    path=str(ckpt_dir / "best_model.pt"), map_location=device)
                if ckpt is not None:
                    unwrapped = model.module if hasattr(model, 'module') else model
                    unwrapped.load_state_dict(ckpt['model_state_dict'])
                    for pg in optimizer.param_groups:
                        pg['lr'] *= 0.5
                    if is_main:
                        logger.info(f"Recovered from best checkpoint, "
                                    f"new lr={optimizer.param_groups[0]['lr']:.2e}")
                consecutive_nan_epochs = 0
            continue
        else:
            consecutive_nan_epochs = 0

        # Validate
        val_loss, per_output = validate(model, val_loader, clinical_weights, device)

        # Scheduler step
        scheduler.step()
        lr_now = optimizer.param_groups[0]['lr']

        # Aggregate across ranks
        if world_size > 1:
            metrics = torch.tensor([train_loss, val_loss], device=device)
            dist.all_reduce(metrics, op=dist.ReduceOp.AVG)
            train_loss, val_loss = metrics[0].item(), metrics[1].item()

        # Record history
        current_window = train_ds.window_size if args.curriculum else args.window_size
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['lr'].append(lr_now)
        history['window_size'].append(current_window)

        elapsed = time.time() - epoch_t0
        is_best = val_loss < best_val_loss
        if is_best:
            best_val_loss = val_loss

        # Logging
        if is_main and (epoch % args.log_every == 0 or epoch == start_epoch or is_best):
            logger.info(
                f"Epoch {epoch:4d}/{args.epochs} | "
                f"train={train_loss:.5f} val={val_loss:.5f} | "
                f"lr={lr_now:.2e} | window={current_window} | "
                f"{elapsed:.1f}s" + (" *BEST*" if is_best else "")
            )

        # TensorBoard
        if tb_writer is not None:
            tb_writer.add_scalar('Loss/train', train_loss, epoch)
            tb_writer.add_scalar('Loss/val', val_loss, epoch)
            tb_writer.add_scalar('LR', lr_now, epoch)
            tb_writer.add_scalar('Window', current_window, epoch)
            for i, key in enumerate(CLINICAL_KEYS):
                tb_writer.add_scalar(f'ValMSE/{key}', per_output[i].item(), epoch)

        # Checkpointing (rank 0 only)
        if is_main:
            unwrapped = model.module if hasattr(model, 'module') else model
            state = {
                'model_state_dict': unwrapped.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_val_loss': best_val_loss,
                'latent_dim': args.latent_dim,
                'norm_stats': norm_stats,
                'history': history,
                'config': {
                    'latent_dim': args.latent_dim,
                    'window_size': current_window,
                    'n_clinical': N_CLINICAL,
                    'n_critical': N_CRITICAL,
                    'clinical_keys': CLINICAL_KEYS,
                    'curriculum': args.curriculum,
                    'batch_size': args.batch_size,
                    'lr': args.lr,
                },
            }

            # Save periodic checkpoint
            if epoch % args.save_every == 0:
                ckpt_mgr.save_checkpoint(state, epoch, val_loss, is_best=is_best)
            elif is_best:
                ckpt_mgr.save_checkpoint(state, epoch, val_loss, is_best=True)

        # Preemption check
        if ckpt_mgr.should_preempt:
            if is_main:
                logger.warning(f"Preemption signal received at epoch {epoch}")
                unwrapped = model.module if hasattr(model, 'module') else model
                state = {
                    'model_state_dict': unwrapped.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'scheduler_state_dict': scheduler.state_dict(),
                    'best_val_loss': best_val_loss,
                    'latent_dim': args.latent_dim,
                    'norm_stats': norm_stats,
                    'history': history,
                }
                ckpt_mgr.save_preempt_checkpoint(state, epoch, val_loss)
                logger.info("Preempt checkpoint saved — exiting cleanly")
            cleanup_distributed()
            sys.exit(0)

    # ---------------------------------------------------------------
    # Final save
    # ---------------------------------------------------------------
    if is_main:
        unwrapped = model.module if hasattr(model, 'module') else model
        final_state = {
            'model_state_dict': unwrapped.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'val_loss': val_loss,
            'best_val_loss': best_val_loss,
            'latent_dim': args.latent_dim,
            'norm_stats': norm_stats,
            'history': history,
            'config': {
                'latent_dim': args.latent_dim,
                'window_size': current_window,
                'n_clinical': N_CLINICAL,
                'n_critical': N_CRITICAL,
                'clinical_keys': CLINICAL_KEYS,
            },
        }
        ckpt_mgr.save_checkpoint(final_state, args.epochs, val_loss,
                                  is_best=(val_loss <= best_val_loss), tag='final')

        # Save history JSON
        hist_path = run_dir / "training_history.json"
        with open(hist_path, 'w') as f:
            json.dump(history, f, indent=2)

        logger.info(f"\nTraining complete. Best val loss: {best_val_loss:.6f}")
        logger.info(f"Checkpoints: {ckpt_dir}")

        if tb_writer:
            tb_writer.close()

    cleanup_distributed()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="HPC training for Neural ODE cardiorenal surrogate")

    # Data
    parser.add_argument('--data', type=str, required=True,
                        help='Path to HDF5 training data')
    parser.add_argument('--train_frac', type=float, default=0.85,
                        help='Fraction of data for training (default: 0.85)')

    # Model
    parser.add_argument('--latent_dim', type=int, default=LATENT_DIM,
                        help=f'Latent ODE dimension (default: {LATENT_DIM})')

    # Training
    parser.add_argument('--epochs', type=int, default=300,
                        help='Total training epochs (default: 300)')
    parser.add_argument('--batch_size', type=int, default=64,
                        help='Per-GPU batch size (default: 64)')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Peak learning rate (default: 1e-3)')
    parser.add_argument('--weight_decay', type=float, default=1e-5,
                        help='AdamW weight decay (default: 1e-5)')
    parser.add_argument('--warmup_epochs', type=int, default=10,
                        help='LR warmup epochs (default: 10)')
    parser.add_argument('--grad_accum', type=int, default=1,
                        help='Gradient accumulation steps (default: 1)')
    parser.add_argument('--max_grad_norm', type=float, default=10.0,
                        help='Max gradient norm for clipping (default: 10.0)')

    # Windows
    parser.add_argument('--window_size', type=int, default=24,
                        help='Trajectory window size (no curriculum, default: 24)')
    parser.add_argument('--curriculum', action='store_true',
                        help='Enable curriculum learning (ramp window size)')
    parser.add_argument('--window_min', type=int, default=8,
                        help='Min window for curriculum (default: 8)')
    parser.add_argument('--window_max', type=int, default=48,
                        help='Max window for curriculum (default: 48)')

    # Checkpoint / Resume
    parser.add_argument('--resume', type=str, default=None,
                        help='Resume path or "auto" for latest checkpoint')
    parser.add_argument('--run_dir', type=str, default=None,
                        help='Run output directory (default: neural_surrogate/runs/<timestamp>)')
    parser.add_argument('--save_every', type=int, default=10,
                        help='Save checkpoint every N epochs (default: 10)')
    parser.add_argument('--keep_checkpoints', type=int, default=3,
                        help='Keep N most recent periodic checkpoints (default: 3)')

    # Logging
    parser.add_argument('--log_every', type=int, default=5,
                        help='Log every N epochs (default: 5)')

    # DataLoader
    parser.add_argument('--num_workers', type=int, default=2,
                        help='DataLoader workers per process (default: 2)')

    args = parser.parse_args()

    # Default run directory
    if args.run_dir is None:
        # Use SLURM_JOB_ID if available for deterministic naming
        job_id = os.environ.get('SLURM_JOB_ID')
        if job_id:
            run_name = f"slurm_{job_id}"
        else:
            from datetime import datetime
            run_name = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.run_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'runs', run_name)

    return args


if __name__ == '__main__':
    args = parse_args()
    train_hpc(args)
