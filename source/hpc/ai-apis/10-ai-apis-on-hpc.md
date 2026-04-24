---
tags: ai api openai hpc slurm rate-limits security
---

# AI APIs on the Tufts HPC Cluster

**Author:** Kyle Monahan, Research Technology, TTS

---

Researchers often need to process large collections of documents, audio files, or datasets using AI APIs such as OpenAI, Anthropic, or similar services. This guide answers the most common questions about running API-based batch workflows on the Tufts HPC cluster safely and efficiently.

---

## Can I Run API-Based Workflows on the Cluster?

Yes. Tufts HPC compute nodes can make outbound requests to external APIs. Follow these practices to keep your workflow well-behaved on shared infrastructure:

1. **Submit through Slurm** - do not run large batch jobs on the login node. See the [HPC documentation](https://go.tufts.edu/cluster) for how to submit jobs.
1. **Request only the resources your job needs** - if the bulk of compute is happening on the API provider's servers, your job likely needs only modest CPU and memory (e.g., 1 CPU, 4 GB RAM).
1. **Test on a small subset first** - run 5–10 files before launching on thousands. Check timing, cost, and error rates before scaling up.
1. **Respect rate limits** - see the [rate limits](#rate-limits) section below.
1. **Set a spending limit** on your API account before starting any batch job.

---

## Storing Your API Key Securely

Never put your API key as plain text in a job script or Python source file. If the script ends up in version control or shared storage, the key is exposed.

There are four recommended approaches.

### Approach 1: Environment Variable (Simplest)

Export the key in the shell session or script that submits your job:

```bash
export OPENAI_API_KEY="sk-..."
sbatch my_job.sh
```

Slurm inherits environment variables from the submitting shell, so the key is available inside the job without appearing in the script file itself.

### Approach 2: Credentials File (More Durable)

**Step 1** - Create a credentials file once, with restricted permissions:

```bash
echo 'export OPENAI_API_KEY="sk-..."' > ~/.openai_cred
chmod 600 ~/.openai_cred
```

**Step 2** - Source it at the top of your Slurm job script:

```bash
#!/bin/bash
#SBATCH --job-name=openai_batch
#SBATCH --array=0-5999%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

source ~/.openai_cred

python process_file.py --index $SLURM_ARRAY_TASK_ID
```

**Step 3** - Read the key from the environment in your Python script:

```python
import os

api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise EnvironmentError("OPENAI_API_KEY not set. Did you source ~/.openai_cred?")
```

This way your key never appears in Slurm script arguments, `sbatch` command history, or stdout/stderr logs.

### Approach 3: GPG-Encrypted Credentials File

GPG encryption adds a layer of protection beyond file permissions alone. Even if someone gains read access to your home directory, the key is unreadable without your GPG passphrase.

**Step 1** - Check that you have a GPG key, or generate one:

```bash
gpg --list-keys
# If no key exists:
gpg --full-generate-key
```

**Step 2** - Create a plain-text credentials file, encrypt it, then delete the plain-text original:

```bash
echo 'export OPENAI_API_KEY="sk-..."' > ~/.openai_cred
gpg --encrypt --recipient your@tufts.edu ~/.openai_cred
rm ~/.openai_cred        # remove the unencrypted copy
```

This produces `~/.openai_cred.gpg`.

**Step 3** - In your Slurm job script, decrypt the file inline:

```bash
#!/bin/bash
#SBATCH --job-name=openai_batch
#SBATCH --array=0-5999%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

eval "$(gpg --quiet --decrypt ~/.openai_cred.gpg)"

python process_file.py --index $SLURM_ARRAY_TASK_ID
```

The key is decrypted into the job's environment at runtime and never written to disk in plain text.

> **Note:** Non-interactive Slurm jobs need access to your GPG key without a prompt. Run `gpg-agent` in your session before submitting, or configure it to cache your passphrase: `gpg-agent --daemon`. The agent caches your passphrase for a configurable duration so batch jobs can decrypt without interruption.

**Step 4** - Read the key from the environment in your Python script (same as Approach 2):

```python
import os

api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise EnvironmentError("OPENAI_API_KEY not set.")
```

### Approach 4: `.env` File with `python-dotenv`

This approach keeps the secret co-located with your project and out of your Slurm script entirely. It is well suited to Python-only workflows where you control the entry point.

**Step 1** - Create a `.env` file in your project directory with restricted permissions:

```bash
echo 'OPENAI_API_KEY="sk-..."' > ~/my_project/.env
chmod 600 ~/my_project/.env
```

**Step 2** - Add `.env` to your `.gitignore` so it is never committed:

```bash
echo '.env' >> ~/my_project/.gitignore
```

**Step 3** - Install `python-dotenv` in your conda environment:

```bash
pip install python-dotenv
```

**Step 4** - Load the `.env` file at the top of your Python script:

```python
import os
from dotenv import load_dotenv

load_dotenv()  # reads .env from the current working directory

api_key = os.environ.get("OPENAI_API_KEY")
if not api_key:
    raise EnvironmentError("OPENAI_API_KEY not set. Is .env present?")
```

Your Slurm script does not need to reference the key at all:

```bash
#!/bin/bash
#SBATCH --job-name=openai_batch
#SBATCH --array=0-5999%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

cd ~/my_project
python process_file.py --index $SLURM_ARRAY_TASK_ID
```

> **Note:** `load_dotenv()` does not overwrite environment variables that are already set. If you have `OPENAI_API_KEY` exported in your shell, that value takes precedence over what is in `.env`.

---

## Rate Limits

Every API provider enforces rate limits - caps on how many requests or tokens you can send per minute. Exceeding them results in throttling errors (`429 Too Many Requests`).

For OpenAI, rate limits are dynamic and depend on your usage tier. The API includes these headers in every response:

| Header                           | Meaning                                  |
| -------------------------------- | ---------------------------------------- |
| `x-ratelimit-remaining-requests` | Requests remaining in the current window |
| `x-ratelimit-remaining-tokens`   | Tokens remaining in the current window   |

Save these headers alongside your responses, especially when debugging. A simple backoff strategy handles temporary rate limit errors gracefully:

```python
import time
import openai


def call_with_backoff(client, **kwargs):
    """Call the OpenAI API with exponential backoff on rate limit errors."""
    for attempt in range(5):
        try:
            return client.chat.completions.create(**kwargs)
        except openai.RateLimitError:
            wait = 2**attempt
            print(f"Rate limited. Waiting {wait}s before retry {attempt + 1}/5...")
            time.sleep(wait)
    raise RuntimeError("Exceeded maximum retries due to rate limiting.")
```

---

## Sizing Your Slurm Array Job

For standard (non-batch) API endpoints, limit concurrent requests to avoid hitting rate limits. A reasonable starting point is `%20` (20 simultaneous tasks) in your `--array` directive:

```bash
#SBATCH --array=0-5999%20
```

Adjust upward cautiously after confirming your tier supports it. If you are consistently hitting rate limits, reduce concurrency or switch to the provider's batch endpoint (see below).

---

## Using the OpenAI Batch API

For large collections of requests, OpenAI's [Batch API](https://platform.openai.com/docs/guides/batch) is worth considering. It accepts all requests in a single API call rather than thousands of separate calls, and offers higher rate limits and lower cost per token. The tradeoff is that results are returned asynchronously (typically within 24 hours) rather than immediately.

Use the standard API if you need results quickly or interactively. Use the Batch API if throughput and cost matter more than latency.

---

## Checkpointing and Progress Tracking

For jobs processing thousands of files, save each result as it completes rather than collecting everything in memory. This lets you resume from where you left off if a job fails or times out.

```python
import os
import json


def save_result(output_dir, file_index, result, metadata=None):
    """Save one API result to disk immediately after it is received."""
    os.makedirs(output_dir, exist_ok=True)
    record = {"result": result}
    if metadata:
        record["metadata"] = metadata
    path = os.path.join(output_dir, f"{file_index}.json")
    with open(path, "w") as f:
        json.dump(record, f)
```

Saving the API response metadata (model version, token counts, timestamps) alongside the result improves research reproducibility and helps with debugging.

To restart an interrupted job, check which output files already exist before making API calls:

```python
def already_processed(output_dir, file_index):
    return os.path.exists(os.path.join(output_dir, f"{file_index}.json"))
```

---

## Complete Example: Processing a File Collection

The example below combines all of the above into a minimal end-to-end script.

### Slurm Job Script (`run_api_job.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=api_batch
#SBATCH --array=0-5999%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err

source ~/.openai_cred

python process_file.py --index $SLURM_ARRAY_TASK_ID
```

Submit with:

```bash
mkdir -p logs
sbatch run_api_job.sh
```

### Python Script (`process_file.py`)

```python
import argparse
import json
import os
import time

import openai


OUTPUT_DIR = "results"
FILES = sorted(os.listdir("data/"))  # replace with your file list


def already_processed(index):
    return os.path.exists(os.path.join(OUTPUT_DIR, f"{index}.json"))


def call_with_backoff(client, **kwargs):
    for attempt in range(5):
        try:
            return client.chat.completions.create(**kwargs)
        except openai.RateLimitError:
            wait = 2**attempt
            print(f"Rate limited. Waiting {wait}s (attempt {attempt + 1}/5)...")
            time.sleep(wait)
    raise RuntimeError("Exceeded maximum retries.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()

    if already_processed(args.index):
        print(f"Index {args.index} already processed - skipping.")
        return

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise EnvironmentError("OPENAI_API_KEY not set.")

    client = openai.OpenAI(api_key=api_key)

    file_path = os.path.join("data", FILES[args.index])
    with open(file_path) as f:
        content = f.read()

    response = call_with_backoff(
        client,
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": content},
        ],
    )

    result_text = response.choices[0].message.content
    metadata = {
        "model": response.model,
        "usage": response.usage.model_dump(),
        "file": FILES[args.index],
    }

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    with open(os.path.join(OUTPUT_DIR, f"{args.index}.json"), "w") as f:
        json.dump({"result": result_text, "metadata": metadata}, f)

    print(f"Saved result for index {args.index}.")


if __name__ == "__main__":
    main()
```

> **Note:** Replace the `FILES` list, `data/` directory, model name, and prompt with your own. The checkpointing, backoff, and key-handling logic does not need to change.

---

## Summary of Best Practices

| Practice                                                        | Why it matters                                 |
| --------------------------------------------------------------- | ---------------------------------------------- |
| Submit through Slurm                                            | Protects shared login nodes                    |
| Request minimal resources                                       | API calls are lightweight locally              |
| Test on 5–10 files first                                        | Catch bugs and estimate cost before full run   |
| Store keys securely (env var, credentials file, GPG, or `.env`) | Prevents accidental key exposure               |
| Save results incrementally                                      | Resume without reprocessing completed files    |
| Save response metadata                                          | Supports reproducibility and debugging         |
| Add a backoff function                                          | Handles transient rate limit errors gracefully |
| Cap concurrent array tasks                                      | Stays within API rate limits                   |
