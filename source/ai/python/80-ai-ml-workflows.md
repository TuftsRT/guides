---
tags: python hpc machine-learning pytorch tensorflow
---

# AI and Machine Learning Workflows on the Tufts HPC

This guide covers how to set up and run deep learning and AI workflows on the Tufts HPC cluster, including PyTorch, TensorFlow, and Hugging Face Transformers. It assumes you are comfortable with conda environments (see [Package Management](60-package-management.md)) and basic SLURM job submission (see the [SLURM documentation](../../hpc/slurm/index.md)).

---

## GPU Access on the Tufts HPC

Deep learning frameworks are much faster on GPUs than CPUs for most training tasks. The Tufts HPC has GPU nodes available on the `gpu` partition.

To request a GPU in an interactive session:

```bash
srun -p gpu --gres=gpu:1 --mem=16g -n 4 -t 2:00:00 --pty bash
```

To request a GPU in a SLURM batch script, add these lines to your `#SBATCH` directives:

```bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
```

Check which GPU types are available:

```bash
sinfo -p gpu -o "%N %G"
```

---

## Installing PyTorch

PyTorch is the most widely used deep learning framework in academic research. Install it with conda, using the official PyTorch install command for your CUDA version.

First, check the CUDA version available on the HPC GPU nodes:

```bash
nvidia-smi
```

Then install PyTorch with the matching CUDA version. For example, for CUDA 12.1 (substitute the version reported by `nvidia-smi` on your node):

```bash
conda create --name pytorch-env python=3.11
conda activate pytorch-env
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
```

> **Always use the official PyTorch installation command** from [pytorch.org/get-started](https://pytorch.org/get-started/locally/). The exact command depends on your CUDA version and may change with new releases. The `pytorch-cuda=12.1` value in the example above is a placeholder; substitute the version reported by `nvidia-smi` on a GPU node.

Verify that PyTorch can see the GPU:

```python
import torch

print(torch.cuda.is_available())  # Should print True on a GPU node
print(torch.cuda.get_device_name(0))
```

---

## Installing TensorFlow

```bash
conda create --name tf-env python=3.11
conda activate tf-env
module load cuda/12.9.0
pip install tensorflow
```

Verify GPU access:

```python
import tensorflow as tf

print(tf.config.list_physical_devices("GPU"))
```

> **CUDA module required:** TensorFlow detects the GPU at runtime via the CUDA libraries. You must load the CUDA module (`module load cuda/12.9.0`) before running TensorFlow, both during installation and in any batch script that uses it. If you omit this, `tf.config.list_physical_devices('GPU')` will return an empty list. If you encounter other GPU errors, check the [TensorFlow GPU guide](https://www.tensorflow.org/install/gpu) for version compatibility.

---

## Hugging Face Transformers

[Hugging Face Transformers](https://huggingface.co/docs/transformers) provides pre-trained models for NLP (BERT, GPT-2, LLaMA, etc.), computer vision, and more.

### Installation

```bash
conda activate pytorch-env
pip install transformers datasets accelerate
```

### Running Inference with a Pre-Trained Model

```python
from transformers import pipeline

classifier = pipeline("sentiment-analysis")
result = classifier("The results were better than expected.")
print(result)
```

> **Model downloads on the HPC:** Hugging Face models are downloaded from the internet on first use. Compute nodes may have restricted internet access, so download models during an interactive session (e.g., via OnDemand shell or a `srun` session) before submitting batch jobs. By default, models are cached in `~/.cache/huggingface/`. If your home directory quota is limited, redirect the cache to lab storage by setting the environment variable:
>
> ```bash
> export HF_HOME=/cluster/tufts/XXXXlab/$USER/huggingface
> ```
>
> Add this line to your batch scripts and your shell profile to keep it consistent.

### Fine-Tuning a Model

```python
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    Trainer,
    TrainingArguments,
)
from datasets import load_dataset

# Load tokenizer and model
model_name = "distilbert-base-uncased"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)

# Load and tokenize dataset
dataset = load_dataset("imdb")


def tokenize(batch):
    return tokenizer(batch["text"], truncation=True, padding=True)


tokenized = dataset.map(tokenize, batched=True)

# Training arguments
training_args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=16,
    eval_strategy="epoch",
    fp16=True,  # Use mixed precision on GPU
)

# Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized["train"],
    eval_dataset=tokenized["test"],
)
trainer.train()
```

---

## SLURM Batch Scripts for Deep Learning

For training jobs that will run for hours, submit a batch job rather than using an interactive session.

### Example: PyTorch Training Job

Create a script `train.sh`:

```bash
#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH -t 12:00:00
#SBATCH -o train_%j.log
#SBATCH -J my_training_job

module purge
module load miniforge/25.3.0
source activate pytorch-env

python train.py --epochs 50 --batch-size 64 --lr 1e-4
```

Submit:

```bash
sbatch train.sh
```

Monitor the job:

```bash
squeue --me          # View your queued and running jobs
tail -f train_*.log  # Follow the log output in real time
```

### Example: Hugging Face Training with `accelerate`

For multi-GPU or distributed training, use Hugging Face `accelerate`:

```bash
#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:2
#SBATCH -n 8
#SBATCH --mem=64g
#SBATCH -t 24:00:00
#SBATCH -o finetune_%j.log

module purge
module load miniforge/25.3.0
source activate pytorch-env

accelerate launch --num_processes 2 fine_tune.py
```

---

## Tips for Efficient Training on the HPC

**Use mixed precision training (fp16/bf16):** Halves memory usage and speeds up training on modern GPUs with minimal accuracy loss. In PyTorch:

```python
from torch.cuda.amp import autocast, GradScaler

scaler = GradScaler()
with autocast():
    output = model(input)
    loss = criterion(output, target)
scaler.scale(loss).backward()
scaler.step(optimizer)
scaler.update()
```

**Monitor GPU utilization:** Check that the GPU is being used efficiently. Low GPU utilization often means the data loading is a bottleneck; use `num_workers` in your `DataLoader`:

```python
from torch.utils.data import DataLoader

loader = DataLoader(dataset, batch_size=64, num_workers=4, pin_memory=True)
```

**Save checkpoints:** Long training runs can be interrupted. Save checkpoints regularly:

```python
torch.save(
    {
        "epoch": epoch,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "loss": loss,
    },
    f"checkpoint_epoch_{epoch}.pt",
)
```

**Use the scratch filesystem for large datasets:** The HPC scratch filesystem (`/cluster/scratch/your_username`) provides faster I/O for large datasets than your home directory.

---

## Purpose-Built AI Tools on the HPC

The Tufts HPC also hosts purpose-built AI tools for specific research tasks, including secure LLM chatbot access, speech recognition, and OCR. See [HPC AI Tools](../ai-tools/index.md) for details.

---

## Getting Help

For AI/ML workflows on the HPC, contact **datalab-support@elist.tufts.edu**.
