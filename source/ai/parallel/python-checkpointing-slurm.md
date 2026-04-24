---
tags: python pytorch tensorflow hpc slurm checkpointing machine-learning
---

# Checkpointing and Job Chaining for Long Training Runs

**Author:** Research Technology, TTS

---

## Why Checkpointing Matters on HPC

Training a large model can take hours or days — but SLURM jobs don't run forever. All partitions on the Tufts HPC cluster (`batch`, `gpu`, and `preempt`) enforce a hard **2-day wall-time limit**. Jobs on the `preempt` partition carry an additional risk: they can be interrupted mid-run without warning when a higher-priority job needs the same resources. Without checkpointing, any training progress made up to that point is lost.

A **checkpoint** is a snapshot of your model's state saved to disk at regular intervals during training. It typically includes:

- **Model weights** — the learned parameters
- **Optimizer state** — momentum and adaptive learning rate terms (critical for optimizers like Adam)
- **Training progress** — the current epoch or step number

When your job is interrupted or hits the time limit, the next job picks up from the last checkpoint rather than starting over.

This guide covers checkpointing in PyTorch and TensorFlow, integrating checkpoints with SLURM batch and array jobs, and two patterns for automatically resubmitting jobs until training completes.

---

## Learning Objectives

By the end of this guide, you will be able to:

1. Save and load training checkpoints in PyTorch.
1. Save and load training checkpoints in TensorFlow.
1. Write SLURM batch scripts that resume from a checkpoint, for both single jobs and array jobs.
1. Automatically resubmit jobs using a self-resubmitting script or a dependency-chained launcher.

---

## Checkpointing in PyTorch

PyTorch does not save checkpoints automatically — you write the save and load logic yourself. This gives you full control over what gets saved and when.

### What to Save

A complete checkpoint should include:

| Field                    | Why it matters                                                                                                            |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------- |
| `model.state_dict()`     | The learned weights — the core of your checkpoint                                                                         |
| `optimizer.state_dict()` | Momentum and adaptive learning rate terms; omitting this causes the optimizer to restart cold, which can degrade training |
| `epoch`                  | Where to resume the training loop                                                                                         |
| `best_val_loss`          | Lets you continue tracking the best model across jobs                                                                     |

### Saving a Checkpoint

Save checkpoints to your lab's research storage (`/cluster/tufts/XXXXlab/$USER/`) rather than your home directory. Home directories have a strict quota; checkpoints for large models can easily fill them.

```python
import torch
import os


def save_checkpoint(state, checkpoint_dir, filename="checkpoint.pt"):
    """Save a training checkpoint to disk."""
    os.makedirs(checkpoint_dir, exist_ok=True)
    path = os.path.join(checkpoint_dir, filename)
    torch.save(state, path)
    print(f"Checkpoint saved to {path}")
```

Call this at the end of each epoch (or every N epochs) and whenever a new best validation loss is reached:

```python
# Save every 5 epochs
if (epoch + 1) % 5 == 0:
    save_checkpoint(
        {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "best_val_loss": best_val_loss,
        },
        checkpoint_dir=args.checkpoint_dir,
    )

# Also save the best model separately
if val_loss < best_val_loss:
    best_val_loss = val_loss
    save_checkpoint(
        {
            "epoch": epoch,
            "model_state_dict": model.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
            "best_val_loss": best_val_loss,
        },
        checkpoint_dir=args.checkpoint_dir,
        filename="best_checkpoint.pt",
    )
```

### Loading a Checkpoint

```python
def load_checkpoint(checkpoint_dir, model, optimizer, filename="checkpoint.pt"):
    """Resume from a checkpoint if one exists, otherwise start from scratch."""
    path = os.path.join(checkpoint_dir, filename)
    if not os.path.exists(path):
        print("No checkpoint found — starting from scratch.")
        return 0, float("inf")

    ckpt = torch.load(path, map_location="cpu")
    model.load_state_dict(ckpt["model_state_dict"])
    optimizer.load_state_dict(ckpt["optimizer_state_dict"])
    start_epoch = ckpt["epoch"] + 1
    best_val_loss = ckpt["best_val_loss"]
    print(f"Resumed from epoch {start_epoch}")
    return start_epoch, best_val_loss
```

`map_location="cpu"` ensures the checkpoint loads correctly even if it was saved on a different GPU. PyTorch will move tensors to the correct device when you call `model.to(device)`.

---

### Complete Resumable Training Script

The script below puts `save_checkpoint` and `load_checkpoint` together into a full training loop. Save this as `train.py`. It accepts `--checkpoint-dir` and `--epochs` from the command line so that SLURM scripts can control them without editing the Python file.

```python
import argparse
import os

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------


def save_checkpoint(state, checkpoint_dir, filename="checkpoint.pt"):
    os.makedirs(checkpoint_dir, exist_ok=True)
    path = os.path.join(checkpoint_dir, filename)
    torch.save(state, path)
    print(f"Checkpoint saved: {path}")


def load_checkpoint(checkpoint_dir, model, optimizer, filename="checkpoint.pt"):
    path = os.path.join(checkpoint_dir, filename)
    if not os.path.exists(path):
        print("No checkpoint found — starting from scratch.")
        return 0, float("inf")
    ckpt = torch.load(path, map_location="cpu")
    model.load_state_dict(ckpt["model_state_dict"])
    optimizer.load_state_dict(ckpt["optimizer_state_dict"])
    print(f"Resumed from epoch {ckpt['epoch'] + 1}")
    return ckpt["epoch"] + 1, ckpt["best_val_loss"]


# ---------------------------------------------------------------------------
# Model (replace with your own)
# ---------------------------------------------------------------------------


class SimpleModel(nn.Module):
    def __init__(self, input_dim=64, num_classes=10):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Linear(128, num_classes),
        )

    def forward(self, x):
        return self.net(x)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--checkpoint-dir",
        required=True,
        help="Path to save/load checkpoints (use lab storage)",
    )
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument(
        "--save-every", type=int, default=5, help="Save a checkpoint every N epochs"
    )
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # --- data (replace with your own DataLoaders) ---
    import numpy as np

    rng = np.random.default_rng(42)
    X = torch.tensor(rng.standard_normal((2000, 64)), dtype=torch.float32)
    y = torch.tensor(rng.integers(0, 10, 2000), dtype=torch.long)
    train_loader = DataLoader(
        TensorDataset(X[:1600], y[:1600]), batch_size=64, shuffle=True
    )
    val_loader = DataLoader(TensorDataset(X[1600:], y[1600:]), batch_size=64)

    # --- model and optimizer ---
    model = SimpleModel().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.CrossEntropyLoss()

    # --- resume if checkpoint exists ---
    start_epoch, best_val_loss = load_checkpoint(args.checkpoint_dir, model, optimizer)

    # --- training loop ---
    for epoch in range(start_epoch, args.epochs):
        model.train()
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            loss = criterion(model(X_batch), y_batch)
            loss.backward()
            optimizer.step()

        # validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                val_loss += criterion(model(X_batch), y_batch).item()
        val_loss /= len(val_loader)
        print(f"Epoch {epoch + 1}/{args.epochs} | val_loss: {val_loss:.4f}")

        # save periodically
        if (epoch + 1) % args.save_every == 0:
            save_checkpoint(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "best_val_loss": best_val_loss,
                },
                args.checkpoint_dir,
            )

        # save best
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            save_checkpoint(
                {
                    "epoch": epoch,
                    "model_state_dict": model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "best_val_loss": best_val_loss,
                },
                args.checkpoint_dir,
                filename="best_checkpoint.pt",
            )

    # write a sentinel file so SLURM scripts know training finished
    flag_path = os.path.join(args.checkpoint_dir, "training_complete.flag")
    open(flag_path, "w").close()
    print(f"Training complete. Flag written to {flag_path}")


if __name__ == "__main__":
    main()
```

> **Note:** Replace `SimpleModel`, `train_loader`, and `val_loader` with your own model and data. The checkpoint, resume, and sentinel file logic does not need to change.

---

## Checkpointing in TensorFlow

TensorFlow provides two built-in tools for checkpointing: `tf.train.Checkpoint` with `CheckpointManager` (for custom training loops) and the `ModelCheckpoint` callback (for training with `model.fit()`).

### Option 1: `tf.train.Checkpoint` and `CheckpointManager`

`tf.train.Checkpoint` saves any combination of model, optimizer, and scalar variables. `CheckpointManager` keeps only the most recent N checkpoints on disk to avoid filling your storage quota.

```python
import argparse
import os
import tensorflow as tf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", required=True)
    parser.add_argument("--epochs", type=int, default=50)
    args = parser.parse_args()

    # --- model and optimizer (replace with your own) ---
    model = tf.keras.Sequential(
        [
            tf.keras.layers.Dense(128, activation="relu", input_shape=(64,)),
            tf.keras.layers.Dense(10),
        ]
    )
    optimizer = tf.keras.optimizers.Adam(learning_rate=1e-3)
    loss_fn = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True)

    # track epoch as a tf.Variable so it is saved with the checkpoint
    epoch_var = tf.Variable(0, trainable=False, dtype=tf.int64)

    checkpoint = tf.train.Checkpoint(model=model, optimizer=optimizer, epoch=epoch_var)
    manager = tf.train.CheckpointManager(checkpoint, args.checkpoint_dir, max_to_keep=3)

    # restore the latest checkpoint if one exists
    if manager.latest_checkpoint:
        checkpoint.restore(manager.latest_checkpoint)
        print(f"Resumed from {manager.latest_checkpoint}")
    else:
        print("No checkpoint found — starting from scratch.")

    start_epoch = int(epoch_var)

    # --- training loop ---
    for epoch in range(start_epoch, args.epochs):
        epoch_var.assign(epoch)

        # replace with your actual training step
        with tf.GradientTape() as tape:
            # logits = model(X_batch, training=True)
            # loss = loss_fn(y_batch, logits)
            pass
        # grads = tape.gradient(loss, model.trainable_variables)
        # optimizer.apply_gradients(zip(grads, model.trainable_variables))

        manager.save()
        print(f"Epoch {epoch + 1}/{args.epochs} complete — checkpoint saved.")

    flag_path = os.path.join(args.checkpoint_dir, "training_complete.flag")
    open(flag_path, "w").close()
    print(f"Training complete. Flag written to {flag_path}")


if __name__ == "__main__":
    main()
```

### Option 2: `ModelCheckpoint` Callback

For training with `model.fit()`, the `ModelCheckpoint` callback is simpler. It saves automatically during training; you reload the saved model manually at the start of the next job.

Use this option when you are using a standard Keras training loop. Use `tf.train.Checkpoint` when you have a custom training loop (i.e., you manage gradient tapes yourself).

```python
import argparse
import os
import tensorflow as tf


def build_model():
    model = tf.keras.Sequential(
        [
            tf.keras.layers.Dense(128, activation="relu", input_shape=(64,)),
            tf.keras.layers.Dense(10, activation="softmax"),
        ]
    )
    model.compile(
        optimizer="adam",
        loss="sparse_categorical_crossentropy",
        metrics=["accuracy"],
    )
    return model


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", required=True)
    parser.add_argument("--epochs", type=int, default=50)
    args = parser.parse_args()

    checkpoint_path = os.path.join(args.checkpoint_dir, "checkpoint.keras")

    # resume from checkpoint if it exists, otherwise build fresh
    if os.path.exists(checkpoint_path):
        model = tf.keras.models.load_model(checkpoint_path)
        print(f"Resumed from {checkpoint_path}")
    else:
        print("No checkpoint found — starting from scratch.")
        model = build_model()

    checkpoint_cb = tf.keras.callbacks.ModelCheckpoint(
        filepath=checkpoint_path,
        save_best_only=True,
        monitor="val_loss",
        verbose=1,
    )

    # replace train_dataset and val_dataset with your own data
    # model.fit(
    #     train_dataset,
    #     validation_data=val_dataset,
    #     epochs=args.epochs,
    #     callbacks=[checkpoint_cb],
    # )


if __name__ == "__main__":
    main()
```

> **Note:** `save_best_only=True` means only the model with the best `val_loss` so far is kept. If you want to resume from an arbitrary interruption rather than from the best model, set `save_best_only=False`.

---

## Integrating Checkpoints with SLURM

### Single GPU Batch Job

Save this as `train.sh`. It sets `CHECKPOINT_DIR` to a path in lab storage and passes it to the Python training script. The Python script reads this path to find and save checkpoints.

Replace `XXXXlab` with your actual lab group name and `my_env` with your conda environment name.

```bash
#!/bin/bash -l
#SBATCH -J train_checkpoint             # job name
#SBATCH --time=02-00:00:00              # 2-day wall-time limit (DD-HH:MM:SS)
#SBATCH -p gpu,preempt                  # gpu or preempt, whichever is free first
#SBATCH -N 1
#SBATCH -n 4                            # CPU cores to support the GPU
#SBATCH --mem=64g
#SBATCH --gres=gpu:a100:1              # request 1 A100 GPU
#SBATCH --constraint="a100-80G"        # request the 80 GB variant
#SBATCH --output=train.%j.%N.out
#SBATCH --error=train.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
module load cuda/12.9.0
source activate my_env

# log GPU info for troubleshooting
nvidia-smi

CHECKPOINT_DIR=/cluster/tufts/XXXXlab/$USER/checkpoints/my_experiment

python train.py \
    --checkpoint-dir "$CHECKPOINT_DIR" \
    --epochs 100 \
    --save-every 5

conda deactivate
```

Submit with:

```bash
sbatch train.sh
```

Each time this job is submitted, the training script checks for an existing checkpoint and resumes from it if found.

---

### SLURM Array Jobs

Array jobs run the same script independently for multiple tasks (e.g., one task per cross-validation fold or hyperparameter configuration). Each task must save its checkpoint to a separate directory so they don't overwrite each other.

Save this as `train_array.sh`:

```bash
#!/bin/bash -l
#SBATCH -J train_array
#SBATCH --time=02-00:00:00
#SBATCH -p gpu,preempt
#SBATCH --array=0-4                     # one task per fold (folds 0–4)
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH --gres=gpu:a100:1
#SBATCH --output=train_fold%a.%j.%N.out   # %a = array task index
#SBATCH --error=train_fold%a.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
module load cuda/12.9.0
source activate my_env

nvidia-smi

# Each task saves to its own subdirectory — no collisions between tasks
CHECKPOINT_DIR=/cluster/tufts/XXXXlab/$USER/checkpoints/fold_${SLURM_ARRAY_TASK_ID}

python train_fold.py \
    --checkpoint-dir "$CHECKPOINT_DIR" \
    --fold "${SLURM_ARRAY_TASK_ID}" \
    --epochs 50

conda deactivate
```

Submit with:

```bash
sbatch train_array.sh
```

`${SLURM_ARRAY_TASK_ID}` is automatically set by SLURM to the task index (0, 1, 2, 3, or 4 in this example). Your Python script reads it to determine which fold or configuration to run. If any task is preempted and requeued, it resumes from `checkpoints/fold_<N>/checkpoint.pt` automatically.

In your Python script, read the fold index like this:

```python
parser.add_argument("--fold", type=int, required=True)
# ...
# use args.fold to select the right data split
```

---

## Auto-resubmission

Sometimes you know upfront that training will take longer than the 2-day wall-time limit. Instead of manually resubmitting the job each time, you can automate resubmission using one of two patterns.

Both patterns rely on the checkpoint logic from the PyTorch or TensorFlow sections above. The Python training script writes a `training_complete.flag` sentinel file when it finishes successfully (see `train.py` from the PyTorch section). The SLURM script checks for this file to decide whether to resubmit.

---

### Pattern A: Self-Resubmitting Script

The simplest approach: at the end of the batch script, check whether training is complete. If not, call `sbatch $0` to requeue the same script.

Save this as `train_self_resubmit.sh`:

```bash
#!/bin/bash -l
#SBATCH -J train_checkpoint
#SBATCH --time=02-00:00:00
#SBATCH -p gpu,preempt
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=64g
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint="a100-80G"
#SBATCH --output=train.%j.%N.out
#SBATCH --error=train.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Your_Tufts_Email@tufts.edu

module purge
module load miniforge/25.3.0
module load cuda/12.9.0
source activate my_env

nvidia-smi

CHECKPOINT_DIR=/cluster/tufts/XXXXlab/$USER/checkpoints/my_experiment

python train.py \
    --checkpoint-dir "$CHECKPOINT_DIR" \
    --epochs 200 \
    --save-every 5

# Resubmit only if training is not yet complete
if [ ! -f "$CHECKPOINT_DIR/training_complete.flag" ]; then
    echo "Training not yet complete — resubmitting."
    sbatch "$0"
else
    echo "Training complete — not resubmitting."
fi

conda deactivate
```

Submit once:

```bash
sbatch train_self_resubmit.sh
```

The job will keep resubmitting itself until the Python script writes `training_complete.flag`. If the Python script exits due to an error before writing the flag, the script will still resubmit — add a check on the exit code if you want to stop on failure:

```bash
python train.py --checkpoint-dir "$CHECKPOINT_DIR" --epochs 200 --save-every 5
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "Training script failed with exit code $EXIT_CODE — not resubmitting."
    conda deactivate
    exit $EXIT_CODE
fi

if [ ! -f "$CHECKPOINT_DIR/training_complete.flag" ]; then
    echo "Training not yet complete — resubmitting."
    sbatch "$0"
fi
```

---

### Pattern B: Dependency Chain Launcher

Submit multiple jobs upfront, each one dependent on the previous completing successfully. This makes the full chain visible in `squeue` immediately and is easier to cancel than a self-resubmitting script.

`--dependency=afterok:<jobid>` tells SLURM to start the next job only if the previous one exits with code 0 (success). If any job in the chain fails, the remaining jobs are cancelled automatically.

First, create your standard batch script (e.g., `train.sh` from the SLURM integration section above — no changes needed).

Then save this launcher as `launch_chain.sh`:

```bash
#!/bin/bash
# launch_chain.sh — submit a chain of dependent jobs
# Usage: bash launch_chain.sh <number_of_jobs>
# Example: bash launch_chain.sh 5

NUM_JOBS=${1:-3}

if [ "$NUM_JOBS" -lt 1 ]; then
    echo "Usage: bash launch_chain.sh <number_of_jobs>"
    exit 1
fi

# Submit the first job
JOB_ID=$(sbatch --parsable train.sh)
echo "Submitted job 1: $JOB_ID"
ALL_JOB_IDS="$JOB_ID"

# Chain the remaining jobs
for i in $(seq 2 "$NUM_JOBS"); do
    JOB_ID=$(sbatch --parsable --dependency=afterok:"$JOB_ID" train.sh)
    echo "Submitted job $i: $JOB_ID (starts after previous)"
    ALL_JOB_IDS="$ALL_JOB_IDS $JOB_ID"
done

echo ""
echo "Submitted $NUM_JOBS chained jobs."
echo "To cancel all pending jobs in this chain, run:"
echo "  scancel $ALL_JOB_IDS"
```

Run it:

```bash
bash launch_chain.sh 5   # submit 5 chained jobs (up to 10 days of training)
```

All five jobs appear immediately in `squeue`. Jobs 2–5 show status `PD` (pending) with reason `Dependency` until the previous job completes.

**Choosing between the two patterns:**

|                   | Pattern A (self-resubmitting)                              | Pattern B (dependency chain)                 |
| ----------------- | ---------------------------------------------------------- | -------------------------------------------- |
| Setup             | Modify one script                                          | Two scripts (`train.sh` + `launch_chain.sh`) |
| Visibility        | Only the running job appears in `squeue`                   | All jobs visible immediately                 |
| Cancellation      | `scancel <current_job_id>` (future jobs not yet submitted) | `scancel <all_job_ids>`                      |
| Stops on failure? | Only if you add exit-code check                            | Yes, automatically via `afterok`             |
| Best for          | Open-ended training (unknown number of jobs needed)        | Known upper bound on jobs needed             |

---
