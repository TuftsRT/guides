---
tags: ai generative-ai hpc gpu open-ondemand-ood data-privacy python claude-code
---

# Using the LLM API (Llama Server) on the Cluster

By Uku-Kaspar Uustalu, Peter Nadel, and Kyle Monahan, Research Technology (RT), Tufts Technology Services

The large language model (LLM) application programming interface (API) is an Open OnDemand (OOD) app that launches a private LLM server ([`llama-server`](https://github.com/ggml-org/llama.cpp/tree/master/tools/server) from llama.cpp) on a graphics processing unit (GPU) node of the Tufts High Performance Compute (HPC) Cluster.

The server is compatible with the [OpenAI API](https://platform.openai.com/docs/api-reference), so open-weight text and vision models can be queried from scripts, notebooks, and applications with no external API keys or services. For a browser chat interface instead, see the [Research Chatbot](./ollama-owui.md).

```{attention}
An HPC account is required to access this application. Please visit the [Account Request](../../hpc/access/10-account-request.md) page to request one.
```

## Data Privacy

Inference runs entirely on Tufts-owned hardware, so by default your prompts and outputs are not sent externally. No request content is logged to an external service, and nothing sent to the server is used to train external models unless you specifically ask. Before using regulated or restricted data, confirm with the relevant Institutional Review Board (IRB) and the [Tufts Privacy Office](https://privacy.tufts.edu/) that this workflow meets the data's requirements. See also the [Guidelines for Use of Generative AI Tools](https://it.tufts.edu/guidelines-use-generative-ai-tools).

```{warning}
The server is protected only by a per-session API key and is intended for use **from the Cluster or the Tufts virtual private network (VPN) only**. Never expose it to the public internet or share the API key.
```

## Getting Started

1. Log into [Open OnDemand](https://ondemand-prod.pax.tufts.edu/) and select **Interactive Apps → LLM API**.
1. Fill out the launch form and click "Launch". Most fields can be left at their defaults:
   - _Number of hours_: How long the server stays up (default: one hour).
   - _Resource preset_: The GPU, RAM, and CPU bundle. For example, `a100-40G/94/4` is an A100 GPU with 40 GB of video memory (VRAM), 94 GB of RAM, and four cores. Larger models need more VRAM.
   - _Use router mode_ (recommended): Loads any model from the shared library on demand, so models can be switched between requests without relaunching.
   - _Use random API key_ (recommended): Generates a new key for each session.
1. Wait for the job to move from _Queued_ to _Running_, then click "Connect to LLM API" to view the session's **host**, **port**, and **API key**.

### Choosing a Model

The shared library at `/cluster/tufts/datasets/ai/models/inference/gguf` includes more than 40 models, and the list changes over time. Each model is a folder whose name is the model ID used in requests, so the current list can be seen from any Cluster shell (e.g., in Open OnDemand under **Clusters → Tufts HPC Shell Access**) with:

```bash
ls /cluster/tufts/datasets/ai/models/inference/gguf
```

Once a session is running, the same list is also available from the server itself (see [curl](#curl) below). Good starting points are:

- **General text:** Qwen3.5 (e.g., `Qwen3.5-27B-Q8_0`), GLM, or gpt-oss
- **Code:** Devstral (e.g., `Devstral-Small-2-24B-Instruct-2512-Q8_0`) or Qwen3-Coder-Next
- **Images:** Qwen3-VL, gemma-4, or Qwen3.5

Start with a smaller model to test a workflow, then move to a larger one if needed. To request that a model be added to the shared library, email {{ email }}.

### How Many Models Fit?

In router mode, all loaded models share the GPU's VRAM. Each model needs at least its size on disk, plus memory for its context window (the text it is working on). Each model is a folder in the shared library, and its size can be checked with `du -sh` (e.g., `du -sh /cluster/tufts/datasets/ai/models/inference/gguf/Qwen3.5-27B-Q8_0`), or estimated as the number of parameters times about 1 byte for 8-bit (`Q8_0`) models or about 0.6 bytes for 4-bit (e.g., `Q4_K_XL` or `q4_0`) models:

| Model                             | Parameters | Size on disk |
| --------------------------------- | ---------- | ------------ |
| `Qwen3.5-27B-Q8_0` (8-bit)        | 27 billion | 28 GB        |
| `gemma-4-31B-it-Q8_0` (8-bit)     | 31 billion | 32 GB        |
| `gemma-4-31B-it-qat-q4_0` (4-bit) | 31 billion | 18 GB        |

The server sizes each model's context window to fit the available VRAM, so a single model can use most of a GPU. In testing, `Qwen3.5-27B-Q8_0` (28 GB on disk) used about 39 GB of a 40 GB A100. In practice, plan on **one model at a time** on a 40 GB GPU, and choose a preset with more VRAM (e.g., an H200 with 141 GB) to work with larger or multiple models. A model that does not fit will either fail to load or run much more slowly with part of it on the CPU.

## Connecting to the Server

Replace `<HOST>`, `<PORT>`, and `<API_KEY>` below with the values from the session page. These change with every launch.

### curl

List the available model IDs, then send a chat message:

```bash
curl http://<HOST>:<PORT>/v1/models -H "Authorization: Bearer <API_KEY>"

curl http://<HOST>:<PORT>/v1/chat/completions \
  -H "Authorization: Bearer <API_KEY>" \
  -H "Content-Type: application/json" \
  -d '{"model": "Qwen3.5-27B-Q8_0", "messages": [{"role": "user", "content": "Hello!"}]}'
```

### Python

Point the [OpenAI Python client](https://github.com/openai/openai-python) (`pip install openai`) at the server:

```python
from openai import OpenAI

client = OpenAI(base_url="http://<HOST>:<PORT>/v1", api_key="<API_KEY>")

response = client.chat.completions.create(
    model="Qwen3.5-27B-Q8_0",
    messages=[{"role": "user", "content": "Hello!"}],
)
print(response.choices[0].message.content)
```

Add `stream=True` to receive the response token by token. Vision models (e.g., `Qwen3-VL-32B-Instruct-Q8_0`) accept images as base64 data URIs:

```python
import base64

with open("figure.png", "rb") as f:
    image_data = base64.b64encode(f.read()).decode("utf-8")

response = client.chat.completions.create(
    model="Qwen3-VL-32B-Instruct-Q8_0",
    messages=[
        {
            "role": "user",
            "content": [
                {"type": "text", "text": "What is in this image?"},
                {
                    "type": "image_url",
                    "image_url": {"url": f"data:image/png;base64,{image_data}"},
                },
            ],
        }
    ],
)
```

### Claude Code

[Claude Code][claude-code-setup-url] can use the server through its Anthropic-compatible `/v1/messages` endpoint. Set these environment variables before running `claude`, using a model ID from `/v1/models`:

```bash
export ANTHROPIC_BASE_URL="http://<HOST>:<PORT>"
export ANTHROPIC_AUTH_TOKEN="<API_KEY>"
export ANTHROPIC_MODEL="<MODEL_ID>"
export ANTHROPIC_DEFAULT_HAIKU_MODEL="<MODEL_ID>"
export CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1  # disable telemetry to Anthropic
```

No Anthropic account is needed, and model traffic stays on the Cluster. File edits and commands require a model that supports tool calling; if they fail, try a different model. To revert, `unset` these variables. For using Claude with Anthropic's cloud, see [Securing Claude Team for Use at Tufts](../getting-started/20-securing-claude-team.md).

## Processing Many Files

The script below sends every `.txt` file in a folder to the model and saves each result as soon as it is returned. It reads the server address and key from the `LLM_BASE_URL` and `LLM_API_KEY` environment variables, so the same script works with an Open OnDemand session or a [batch job](#advanced-running-the-server-via-slurm). Save it as `process_files.py` and run it on the Cluster from a Python environment with the `openai` package installed (see [Python for Research and Data Science](../python/index.md)):

```python
import json
import os
from pathlib import Path

from openai import OpenAI

MODEL = "Qwen3.5-27B-Q8_0"
scratch = Path("/cluster/scratch") / os.environ["USER"]
input_dir = scratch / "documents"
output_file = scratch / "results.jsonl"

client = OpenAI(base_url=os.environ["LLM_BASE_URL"], api_key=os.environ["LLM_API_KEY"])

# Skip files already saved by an earlier run
done = set()
if output_file.exists():
    with open(output_file) as f:
        done = {json.loads(line)["file"] for line in f}

for path in sorted(input_dir.glob("*.txt")):
    if path.name in done:
        continue
    response = client.chat.completions.create(
        model=MODEL,
        messages=[{"role": "user", "content": f"Summarize:\n\n{path.read_text()}"}],
        temperature=0,
        seed=42,
    )
    result = {
        "file": path.name,
        "model": MODEL,
        "output": response.choices[0].message.content,
    }
    with open(output_file, "a") as f:
        f.write(json.dumps(result) + "\n")
```

For an Open OnDemand session, set the variables first (`<HOST>`, `<PORT>`, and `<API_KEY>` come from the session page):

```bash
export LLM_BASE_URL="http://<HOST>:<PORT>/v1"
export LLM_API_KEY="<API_KEY>"
python process_files.py
```

The results file can be loaded for analysis with `pandas.read_json("results.jsonl", lines=True)`.

### Speeding Up Large Jobs

- **Turn off thinking for simple tasks.** Reasoning models such as Qwen3.5 "think" before answering by default, which can generate hundreds of hidden tokens per request. For straightforward tasks like summarizing or classifying, turn thinking off by adding `extra_body={"chat_template_kwargs": {"enable_thinking": False}}` to `client.chat.completions.create()`. In testing, this cut a one-sentence summary from 32 seconds to under one second with the same quality of answer.
- **Send requests in parallel.** By default, the server can process four requests at once. Sending several requests at the same time (e.g., with Python's `concurrent.futures.ThreadPoolExecutor`) can speed up large jobs.

### Checkpointing and Saving to Scratch

The script above checkpoints its progress: each result is written to `results.jsonl` as soon as it arrives, and files already in the results are skipped. If the session ends, runs out of time, or is preempted partway through, run the script again in a new session to continue where it left off.

Both the input documents and the results are kept in the scratch directory (`/cluster/scratch/<username>`), which is designed for large, temporary working data. Scratch has a 100 GB quota per user, and **files not accessed in 21 days are automatically deleted**, so copy final results to home or project storage when finished. See [HPC Storage](../../hpc/access/30-storage.md) for details.

### Reproducible Results

For research use, record enough information to reproduce or report the results:

- Set `temperature=0` and a fixed `seed` (as in the script above) so that the same input gives the same output as closely as possible.
- Record the exact model ID, including its quantization (e.g., `Q8_0`), and the `llama.cpp` module version. The script saves the model ID with each result.
- Outputs can still vary slightly between GPU types and software versions, and models can be wrong. Review a sample of the outputs before relying on them.

## Limits and Stopping the Session

The LLM API has no per-token charges; it uses the normal HPC allocation like any GPU job, for as long as the session runs. Because GPUs are shared, delete the session from the Open OnDemand "My Interactive Sessions" page as soon as it is no longer needed. Otherwise, it ends when the requested hours run out. Requests that switch between models trigger reloads, which adds latency.

## Advanced: Running the Server via Slurm

For large jobs, a single [batch job](../../hpc/slurm/batchjob.md) can start the server, run `process_files.py` (see [Processing Many Files](#processing-many-files)), and release the GPU when it finishes, with no Open OnDemand session to watch. Save the following as `llm-batch.sh` and submit it with `sbatch llm-batch.sh`:

```bash
#!/bin/bash -l
#SBATCH -J llm-batch
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --constraint="a100-40G"
#SBATCH --cpus-per-task=4
#SBATCH --mem=94g
#SBATCH --time=04:00:00
#SBATCH --output=llm-batch.%j.out

module purge
module load ai-tools "llama.cpp/0.6.0"

# Random key and port for this job; the server only listens on this node
export LLAMA_API_KEY=$(openssl rand -hex 32)
PORT=$(shuf -i 20000-40000 -n 1)

llama-server --threads 4 --offline --host 127.0.0.1 --port $PORT \
  --models-dir /cluster/tufts/datasets/ai/models/inference/gguf &
SERVER_PID=$!

# Wait until the server is ready (stop if it fails to start)
until curl -sf http://127.0.0.1:$PORT/health > /dev/null; do
  kill -0 $SERVER_PID 2> /dev/null || exit 1
  sleep 5
done

# Python environment with the openai package installed
module load miniforge/25.3.0
source activate <ENV_NAME>

export LLM_BASE_URL="http://127.0.0.1:$PORT/v1"
export LLM_API_KEY="$LLAMA_API_KEY"
python process_files.py

kill $SERVER_PID
```

Adjust the GPU (`--constraint`), memory, and time to fit the model and workload. Because `process_files.py` checkpoints its progress, a job that runs out of time can simply be resubmitted. Run `module load ai-tools` and then `module avail llama.cpp` to check for newer versions.

## Troubleshooting

- **Job stuck in the queue:** Choose a smaller GPU preset or a wider VRAM range.
- **Connection refused:** Confirm the job is _Running_, the host, port, and key match the current session, and the client is on the Cluster or VPN.
- **Slow first response:** Large models can take a few minutes to load.
- **Model not found or out of memory:** Model IDs match the folder names in the shared library; check the ID against `/v1/models`. For memory errors, relaunch with a larger GPU or choose a smaller model.

For any questions, please reach out to Research Technology at {{ email }}.

[claude-code-setup-url]: https://code.claude.com/docs/en/setup
