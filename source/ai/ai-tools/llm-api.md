---
tags: ai hpc gpu open-ondemand-ood
---

# Using the LLM API on the Cluster

By Uku-Kaspar Uustalu, Peter Nadel, and Kyle Monahan, Research Technology (RT), Tufts Technology Services

The LLM API application (listed in Open OnDemand as "LLM Server (LLaMA CPP)") launches an [OpenAI API](https://github.com/openai/openai-openapi)-compatible large language model (LLM) server on a dedicated graphics processing unit (GPU) allocation on the Tufts High Performance Compute (HPC) Cluster. Unlike the [Research Chatbot](./ollama-owui.md), which provides a chat window for open-ended, human-driven conversation, the LLM API is meant to be called _programmatically_: from a coding assistant such as [Continue.dev](https://continue.dev) in VS Code, from a script or notebook using `curl` or the [OpenAI Python client](https://github.com/openai/openai-python), or from any other tool that already knows how to speak to an OpenAI-compatible endpoint.

Because the server runs entirely on Tufts-managed hardware for the life of a single job, none of the prompts, code, or data sent to it leave the Cluster, unless specifically downloaded or shared by the user. (For more information on the secure use of generative AI tools at Tufts, please consult the Tufts Technology Services [Guidelines for Use of Generative AI Tools](https://it.tufts.edu/guidelines-use-generative-ai-tools).)

This document does not provide an in-depth description of every configuration option. It focuses on the choices a new user actually needs to make and on how to connect to a running server once it is up.

## Who is this for?

The LLM API is best suited for researchers and developers who want to plug a large language model into a tool or workflow rather than chat with it directly: for example, using an AI coding assistant inside VS Code, calling an LLM from a Python notebook or pipeline, or prototyping an application against an OpenAI-compatible endpoint before pointing it at a commercial API. If a simple chat interface is all that is needed, the [Research Chatbot](./ollama-owui.md) is likely a better fit.

## What is the LLM API?

The LLM API is Tufts' implementation of the [llama.cpp server](https://github.com/ggml-org/llama.cpp/tree/master/tools/server) on the HPC Cluster. llama.cpp is an open-source program for running LLMs that have been converted to the GGUF (GPT-Generated Unified Format) file format, a compact, quantized format that allows large models to run efficiently on a single GPU. The llama.cpp server exposes an HTTP API that mirrors the OpenAI API, so any tool built to talk to OpenAI's API (`/v1/chat/completions`, `/v1/models`, and so on) can talk to it with only the server address and an API key changed.

For a deeper look at the server itself, see the [llama.cpp server documentation](https://github.com/ggml-org/llama.cpp/tree/master/tools/server). For the API shape it implements, see the [OpenAI API reference](https://github.com/openai/openai-openapi).

## Getting started

The LLM API is an Open OnDemand (OOD) application on the Cluster, meaning it can be accessed from the Interactive Apps drop-down menu in the Open OnDemand website. To get started, visit and log into the [Open OnDemand website for the Tufts Cluster](https://ondemand-prod.pax.tufts.edu/). Once there, select the "Interactive Apps" drop-down menu and click on "LLM Server (LLaMA CPP)".

## Configuring your session

The launch form has a lot of fields, but most users only need to think about a handful of them. Everything under "Advanced" below can safely be left at its default value.

### Job runtime

- _Hours_: How long the session will run before it is automatically terminated, up to a maximum of 48 hours. Jobs requesting more hours may take longer to start. Choose a value with some buffer -- an idle session can simply be ended early if it finishes ahead of schedule.

### Compute resources

- _Resource Preference_: Choose "Use Resource Presets (Recommended)" unless there is a specific reason to hand-pick every resource. Presets bundle a GPU, an amount of video memory (VRAM, the GPU's own memory used to hold the model), system RAM, and CPU cores into a single, tested combination.
  - _Resource Preset_: The available presets range from an A100 GPU with 40 GB VRAM up to an H200 GPU with 140 GB VRAM. Pick a preset whose VRAM can hold the model to be used -- models that do not fully fit into VRAM will run more slowly, and models larger than the combined VRAM and RAM cannot run at all.
  - _Preset Multiplier_: Multiplies the whole preset (GPUs, VRAM, RAM, and CPU cores together) by an integer, for cases where a single preset is not enough (for example, a very large model split across multiple GPUs).
  - "Manually Specify Configuration (Advanced)" instead exposes individual controls for CPU cores, RAM, GPU preference (no GPU, any GPU, a minimum or exact VRAM requirement, or a specific GPU model), and number of GPUs. "Allow Preemption" is also only shown here, and trades a potentially faster job start for the risk that the job could be stopped without warning if higher-priority work needs the same resources.

### Server mode and models

- _Server Mode_: Choose "Router (Multi-Model)" for most use cases. Router mode watches a directory of models and loads or unloads them on demand, so a model can be switched from a client (e.g. Continue.dev) without cancelling and resubmitting the job. "Classic (Single-Model)" preloads exactly one model for the entire lifetime of the job; it cannot be changed without ending the session and starting a new one, but it is the mode required to run a custom GGUF file that is not part of the curated, pre-downloaded set.
- _Model Set_ (Router mode only): "Default (Recommended)" points at a curated, pre-downloaded directory of GGUF models maintained by RT (currently including families such as Devstral, GLM, Kimi, MiniMax, Nemotron, several generations of Qwen, Gemma, gpt-oss, and Granite, among others). "Custom (Advanced)" allows a Models Directory and, optionally, a Model Presets File (`.ini`) to be specified instead -- see the [llama.cpp server model sources](https://github.com/ggml-org/llama.cpp/tree/master/tools/server#model-sources) and [model presets](https://github.com/ggml-org/llama.cpp/tree/master/tools/server#model-presets) documentation for the expected file layout.
- _Model Preference_ (Classic mode only): "Select Pre-Downloaded Model (Recommended)" offers the same curated set as a dropdown. "Manually Specify Custom Model File (Advanced)" instead asks for a path to a GGUF file on the Cluster (for a sharded model, point to the first shard, the one ending in `-00001-of-?????.gguf`, and keep every shard in the same directory) and, optionally, a Multi-Modal Projector (mmproj) file, needed only for models with vision or other multi-modal capabilities and usually distributed alongside the model itself.

The pre-downloaded model catalog changes over time as models are added; do not rely on any specific list from this document. The current set is always visible in the OOD form's model dropdown, and, once a Router-mode session is running, the exact model IDs available to that session can be queried with `GET /v1/models` (see [Connecting to your running session](#connecting-to-your-running-session) below).

### API key

- _API Key Preference_: "Randomly Generate (Recommended)" creates a strong key automatically and displays it once the session starts. "Manually Specify API Key (Advanced)" allows a specific key of at least 32 alphanumeric characters to be supplied instead, for example to keep the same key across several relaunches.

### Advanced server tuning (optional)

The remaining fields configure the underlying llama.cpp server directly and are safe to leave at their defaults unless there is a specific, known reason to change them:

- _Context Size_, _Tokens to Predict_, _Logical/Physical Batch Size_ control the prompt window size, how many tokens are generated per response, and internal batching. Leaving Context Size and Tokens to Predict at `0`/`-1` defers to the model's own defaults.
- _RoPE Parameters_ and the associated _YaRN_ settings are ways of stretching a model's context length beyond what it was originally trained on (RoPE stands for Rotary Position Embedding; YaRN is one extrapolation method for it). These are advanced tuning knobs -- leave them at their defaults unless a specific model calls for extended context and the appropriate values are already known.
- _Enable KV Cache Offloading_ and the _KV Cache Data Type_ options control where the attention key/value cache lives (GPU vs. system memory) and at what numeric precision. Lower-precision types (e.g. `q8_0`, `q4_0`) reduce VRAM use at some cost to output quality.
- _GPU Offload Strategy_ controls how many of a model's layers are placed on the GPU. "Automatic" works for most cases; manual layer counts are mainly useful for squeezing a large model into limited VRAM.

## Connecting to your running session

Once the session reaches the running state, Open OnDemand displays a connection card with everything needed to reach the server.

### Reading the connection card

- _Base URL_: The server's address, in the form `http://<host>.pax.tufts.edu:<port>`. This is the root URL for every API call; the OpenAI-style endpoints live under `<Base URL>/v1/...`.
- _Model ID_: In Classic mode, this is the single fixed model's filename. In Router mode, no single Model ID is shown -- the card notes that multiple models are available and that they can be listed by calling `GET /v1/models` on the running server.
- _API Key_: The key generated (or specified) at launch time. It is required by every request as a bearer token.
- A ready-to-paste Continue.dev configuration block (see below).
- A "Launch Web UI" button, which opens llama.cpp's built-in browser chat interface -- useful for a quick manual test without configuring any external tool.
- A "View Server Log" link to the running job's live log output, useful for troubleshooting (see [Troubleshooting](#troubleshooting) below).

### VS Code with Continue.dev

[Continue.dev](https://continue.dev) is a VS Code (and JetBrains) extension that adds an AI coding assistant backed by a model of choice. Copy the block from the connection card into the `models` list of Continue's `config.yaml` (accessible from the Continue panel in VS Code, or at `~/.continue/config.yaml`):

```yaml
models:
  - name: Tufts LLM Endpoint
    provider: openai
    model: AUTODETECT
    apiBase: http://<host>.pax.tufts.edu:<port>/v1
    apiKey: <your-api-key>
    roles:
      - chat
      - edit
      - apply
```

Replace `<host>`, `<port>`, and `<your-api-key>` with the values from the connection card. In Router mode, `model: AUTODETECT` lets Continue discover whatever model(s) the server is currently serving; in Classic mode, replace `AUTODETECT` with the exact Model ID shown on the card.

### curl

For quick tests or scripting outside of Continue.dev, the server can be called directly with `curl`. To send a chat message:

```bash
curl http://<host>.pax.tufts.edu:<port>/v1/chat/completions \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer <your-api-key>" \
  -d '{
    "model": "<model-id>",
    "messages": [
      { "role": "user", "content": "Hello, how are you?" }
    ]
  }'
```

In Router mode, `<model-id>` should match one of the models returned by:

```bash
curl http://<host>.pax.tufts.edu:<port>/v1/models \
  -H "Authorization: Bearer <your-api-key>"
```

In Classic mode, the server is only ever serving the one fixed model, so the `model` field can typically be left as the Model ID shown on the connection card.

### OpenAI Python client

The official [OpenAI Python client](https://github.com/openai/openai-python) can be pointed at the server by overriding its `base_url` and `api_key`:

```python
from openai import OpenAI

client = OpenAI(
    base_url="http://<host>.pax.tufts.edu:<port>/v1",
    api_key="<your-api-key>",
)

response = client.chat.completions.create(
    model="<model-id>",
    messages=[{"role": "user", "content": "Hello, how are you?"}],
)

print(response.choices[0].message.content)
```

This same pattern works from a Jupyter notebook running on or off the Cluster, as long as it can reach the Base URL.

### Claude Code

[Claude Code](https://claude.com/product/claude-code) is Anthropic's command-line coding agent. Unlike Continue.dev or the OpenAI Python client, Claude Code speaks the [Anthropic Messages API](https://docs.anthropic.com/en/api/messages) rather than the OpenAI API. Recent versions of the llama.cpp server used by the LLM API expose an Anthropic-compatible `/v1/messages` endpoint alongside the OpenAI-compatible routes described above, which allows Claude Code to be pointed directly at a running LLM API session without any third-party translation proxy.

```{attention}
Tool calling (used by Claude Code for file edits and running commands) requires a model whose chat template supports tool use. This is the case for all models currently available on the Tufts HPC, but individual tools may have unexpected results, and this has not been independently verified for every model and tool combination. If Claude Code connects but cannot edit files or run commands, contact {{ email }} to confirm tool-calling support for the session's model.
```

#### Installing Claude Code

Install Claude Code by following the [official installation guide][claude-code-setup-url], which covers the official installer as well as Homebrew, WinGet, npm, and Linux package manager alternatives.

#### Connecting Claude Code to the LLM API

Point Claude Code at the running session instead of Anthropic's cloud by setting the following environment variables before launching `claude`:

```bash
export ANTHROPIC_BASE_URL="http://<host>.pax.tufts.edu:<port>"
export ANTHROPIC_AUTH_TOKEN="<your-api-key>"
export ANTHROPIC_MODEL="<model-id>"
export ANTHROPIC_DEFAULT_HAIKU_MODEL="<model-id>"
```

Replace `<host>`, `<port>`, and `<your-api-key>` with the values from the connection card, and `<model-id>` with an ID returned by `GET /v1/models` (see [curl](#curl) above). `ANTHROPIC_DEFAULT_HAIKU_MODEL` controls the model used for Claude Code's lightweight background tasks; setting it to a model actually hosted by the session avoids errors from Claude Code requesting a default Anthropic model ID the LLM API server does not recognize.

#### Preventing data from reaching Anthropic

Once `ANTHROPIC_BASE_URL` is set, chat and tool-use traffic is sent to the LLM API session rather than to Anthropic, and no Anthropic account or login is required -- `ANTHROPIC_AUTH_TOKEN` can be set to the LLM API session's own API key.

Claude Code independently sends operational telemetry (usage metrics and error reports) to Anthropic and third-party services regardless of where model traffic is routed. To disable this, set the following alongside the connection variables above:

```bash
export CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1
```

```{note}
Anthropic's data-training policy differs by account type: Free, Pro, and Max plans may use conversation data for training unless a user opts out in their [claude.ai privacy settings](https://claude.ai/settings/data-privacy-controls), while Team, Enterprise, and API accounts do not train on data by default. This policy only applies to traffic that reaches Anthropic's servers -- it is not relevant to the LLM API traffic described above, which never leaves Tufts-managed infrastructure. It does apply if Claude Code is later used against Anthropic's cloud (for example, by unsetting `ANTHROPIC_BASE_URL`), so a paid Team, Enterprise, or API account is recommended for that use case.
```

#### Reverting the configuration

To return Claude Code to its default configuration, unset the environment variables set above (or remove them from the shell profile or `~/.claude/settings.json` if they were added there):

```bash
unset ANTHROPIC_BASE_URL ANTHROPIC_AUTH_TOKEN ANTHROPIC_MODEL ANTHROPIC_DEFAULT_HAIKU_MODEL CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC
```

To fully uninstall and reinstall Claude Code instead:

```{warning}
The commands below permanently delete the Claude Code binary and all local configuration, including any settings unrelated to this guide. Double-check each path before running, as `rm -rf` cannot be undone.
```

```bash
rm -f ~/.local/bin/claude
rm -rf ~/.local/share/claude
rm -rf ~/.claude
rm -f ~/.claude.json
```

Reinstalling using the [official installation guide][claude-code-setup-url] restores Claude Code to its default state, prompting for an Anthropic account login on first run (unless the environment variables above are set again).

```{note}
Anyone with the API key shown on the connection card can access the server and read any data sent to it. Do not share it. The key is scoped only to that particular job -- ending the session and launching a new one generates a new key by default (unless a manual key was specified at launch and reused intentionally).
```

## Troubleshooting

- **Where is the server log?** The "View Server Log" link on the connection card opens the running job's `output.log` file, which shows the live llama.cpp server output, including model loading progress and any errors.
- **The server seems unresponsive right after launch.** This is often normal, not a failure. Large GGUF model files are read from network storage, and loading a large quantized model for the first time -- especially immediately after the job starts -- can take a few minutes. Check the server log for loading progress before assuming something is wrong.
- **Requests fail with an authentication error.** Double-check that the API key was copied exactly (it is case-sensitive) and that it matches the key currently shown on the connection card for that job -- keys are not shared between sessions unless manually reused.
- **A Router-mode model cannot be found.** Confirm the model ID matches one of the entries returned by `GET /v1/models`; model IDs are derived from filenames in the models directory and may not exactly match a model's display name.
- **The session ended earlier than expected.** Check whether "Allow Preemption" was enabled in a manually configured job; preemptible jobs can be stopped without warning if higher-priority jobs need the same resources. Otherwise, confirm the requested number of Hours was sufficient.

For any questions, please reach out to Research Technology at: tts-research@tufts.edu.

[claude-code-setup-url]: https://code.claude.com/docs/en/setup.md
