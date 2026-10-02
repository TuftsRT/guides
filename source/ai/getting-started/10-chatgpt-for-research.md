---
tags: ai generative-ai chatgpt codex data-privacy
---

# Using and Installing OpenAI's ChatGPT Locally for Research at Tufts

ChatGPT software, including ChatGPT Work and Codex, can accelerate research and improve productivity. This document provides guidance for Tufts faculty, staff, and students on how to install and use ChatGPT Work (Desktop) and Codex, what data is appropriate for each, and how to protect that data.

For guidance on data classification and when a Tufts Technology Services (TTS) AI review is required, see Using Generative AI and Agentic AI in Research. For how your data moves through OpenAI's systems, see How Generative AI Processes Your Data. See also [Generative AI tools at Tufts][genai-tufts].

## Choosing the Right Tool

Use the ChatGPT Work app for everyday questions, drafting, literature reviews, or working interactively with one prompt at a time. Use ChatGPT Codex, or application programming interface (API) access, when you are automating work such as batch processing, data pipelines, or repeatable programmatic workflows that go beyond a single chat session.

## Which Data to Use

```{warning}
You cannot use confidential or regulated data (Institutional Data Level 3) with ChatGPT Work (Desktop) or Codex.
```

If you're not sure which level your data falls under, check the [Data Finder][data-finder] before you start.

For your research data, you may want to use these products for general writing and coding assistance, literature reviews, or datasets you expect to publish online. For any questions on using ChatGPT products with research data, contact <tts-research@tufts.edu>.

## Using ChatGPT Work (Desktop)

ChatGPT Work is the autonomous agent mode built into the ChatGPT Desktop app. It runs tasks on your machine (or in the cloud, if you choose a remote session) while you do other work, rather than requiring you to sit in a chat window. It is an autonomous agentic tool, and has the potential to read, modify, and share all files or folders on your computer.

### How to Install

- Download the ChatGPT Work app from [openai.com/chatgpt-work](https://openai.com/chatgpt-work/) for Mac or Windows.
- Sign in with your ChatGPT account. Work's autonomous task features require a paid plan (Standard or Pro). A free, Plus, or personal Pro plan allows OpenAI to train on your data by default (via the "Improve the model for everyone" setting) and is not recommended for Tufts research use.
- Once signed in, switch to Work mode to move out of Chat mode and into agentic Work mode. See Protecting Your Data below for how to use it safely.

### Protecting Your Data

- If you pay for a Business plan, OpenAI does not train on your prompts, files, or ChatGPT's replies by default. The exception is if you submit feedback, a bug report, or trigger a sensitive content warning, which sends your current chats and any files within the model context to OpenAI.
- Use Incognito chats for one-off sensitive queries you don't want kept in your regular history.
- Work gives ChatGPT access to a specific folder on your computer. Point it at a scoped working directory rather than your whole file system, and check what it touched when a task finishes. By default, Work will request access to your user folder containing your Documents, Desktop, Downloads, and similar directories. Do not allow this.
- If you connect Work to Google Drive, Slack, or similar, remember that this provides live access to the data in those accounts going forward. You should instead provide limited access to where your public data are stored.

## Using Codex

Codex is the agentic coding tool, available as a terminal command-line interface (CLI), inside the ChatGPT Desktop app's Codex mode, or as an extension in Visual Studio Code and JetBrains.

### How to Install

- Follow the instructions on [OpenAI's website](https://chatgpt.com/codex/) for your operating system of interest.
- After installation, you can confirm the software has installed with `codex --version`.
- Your project directory should only include files you wish to send to Codex.
- You can run `codex` inside a project directory to start your first session.
- Codex requires a paid ChatGPT plan for Tufts research use. A limited free trial exists, but doesn't include the no-training data terms described above.

### Protecting Your Data

Protecting your data in Codex is more technical.

Add explicit deny rules in `config.toml` so Codex can't read sensitive paths, regardless of what you or Codex ask for. To do this, open `config.toml` in the hidden `~/.codex/` directory and add the following TOML:

```toml
approval_policy = "on-request"
default_permissions = "safe-workspace"

# Stop Codex writing session transcripts to disk
history.persistence = "none"

# Disable /feedback submission and usage analytics
feedback.enabled = false
analytics.enabled = false

[permissions.safe-workspace.filesystem]
"/Users/<your-username>/.ssh" = "deny"
"/Users/<your-username>/.aws" = "deny"

[permissions.safe-workspace.filesystem.":workspace_roots"]
"." = "write"
".env" = "deny"
".env.*" = "deny"
"secrets/**" = "deny"
"**/*.pem" = "deny"
```

Do not add `sandbox_mode` to this configuration. OpenAI's reference states that `default_permissions` must not be combined with `sandbox_mode` or `[sandbox_workspace_write]`, because the named permission profile replaces them. Filesystem keys take an absolute path or a special token such as `:workspace_roots`, so spell out your home directory rather than using `~`.

```{important}
Named permission profiles require Codex 0.138.0 or later. Version 0.137.0 and earlier **ignore** `allowed_permission_profiles` and managed `default_permissions` rather than reporting an error, so an older client silently runs without the restrictions above. Check `codex --version` before relying on this configuration.
```

- Session transcripts are saved locally in plaintext. Set `history.persistence = "none"` to stop Codex writing them to `history.jsonl` at all. Server-side retention is governed by your plan's data terms rather than by any client setting.

- `/feedback` sends your conversation history, including code, to OpenAI as a diagnostic report. Codex's documentation does not state whether API keys and tokens are redacted first, so treat anything sent this way as unredacted. Set `feedback.enabled = false` to disable submission across local Codex clients if you're working on anything you wouldn't want leaving your machine.

- Analytics send operational health and usage metrics, not code or file paths. Set `analytics.enabled = false` to turn them off. Raw prompts are never exported with OpenTelemetry logs unless you explicitly opt in via `otel.log_user_prompt`.

- Codex defaults to an approval policy that can auto-approve routine actions without prompting you. For research data we recommend switching to a stricter approval policy by setting `approval_policy = "on-request"` in `config.toml`, so you review each action yourself.

### Protecting Your API Key

If you authenticate Codex with an API key instead of a ChatGPT login, protect that key like a password.

- Never commit keys to a repository. Keep them out of source control. Use a config file like `config.toml`, environment variables, or a secrets manager.
- Never share your key or paste it into tickets, shared documents, or chat.
- Rotate your key periodically, and immediately if you suspect it has been exposed.
- Use one key per person or service where possible, so usage is attributable and a single key can be revoked without disrupting everyone else.
- Scrub keys from logs and screenshots before sharing them for troubleshooting. This matters especially given the session transcript caching and `/feedback` behavior described above.

[data-finder]: https://access.tufts.edu/data-finder
[genai-tufts]: https://it.tufts.edu/ai/generative-ai-tools
