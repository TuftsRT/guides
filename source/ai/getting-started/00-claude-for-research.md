---
tags: ai generative-ai claude claude-code data-privacy
---

# Using and Installing Anthropic's Claude Locally for Research at Tufts

Claude software, including Claude Cowork and Claude Code, can accelerate research and improve productivity. This document provides guidance for Tufts faculty, staff, and students on how to install and use Claude Cowork (Desktop) and Claude Code, what data is appropriate for each, and how to protect that data.

For guidance on data classification and when a Tufts Technology Services (TTS) AI review is required, see Using Generative AI and Agentic AI in Research. For how your data moves through Claude's systems, see How Generative AI Processes Your Data. See also [Generative AI tools at Tufts][genai-tufts].

## Choosing the Right Tool

Use Claude Desktop for everyday questions, drafting, literature reviews, or working interactively with one prompt at a time. Use Claude Code, or application programming interface (API) access, when you are automating work such as batch processing, data pipelines, or repeatable programmatic workflows that go beyond a single chat session.

## Which Data to Use

```{warning}
You cannot use confidential or regulated data (Institutional Data Level 3) with Claude Cowork (Desktop) or Claude Code.
```

If you're not sure which level your data falls under, check the [Data Finder][data-finder] before you start.

For your research data, you may want to use these products for general writing and coding assistance, literature reviews, or datasets you expect to publish online. For any questions on using Claude products with research data, contact <tts-research@tufts.edu>.

## Using Claude Cowork (Desktop)

Claude Cowork is the autonomous agent mode built into the Claude Desktop app. It runs tasks on your machine (or in the cloud, if you choose a remote session) while you do other work, rather than requiring you to sit in a chat window. It is an autonomous agentic tool, and has the potential to read, modify, and share all files or folders on your computer.

### How to Install

- Download Claude Desktop from [claude.ai/download](https://claude.ai/download) for Mac, Windows, or Linux. There's no separate download for Cowork; it's a tab inside the same app.
- Sign in with your Claude account. Cowork's autonomous task features require a paid plan (Pro, Max, Team, or Enterprise). A free plan allows Anthropic to train on your data and is not recommended for Tufts research use.
- Once signed in, click the Cowork tab to switch out of Chat mode and into agentic Cowork mode. See Protecting Your Data below for how to use it safely.

### Protecting Your Data

- If you pay for a Team plan, Anthropic does not train on your prompts, files, or Claude's replies by default. The exception is if you submit feedback, a bug report, or trigger a sensitive content warning, which sends your current chats and any files within the model context to Anthropic.
- Use Incognito chats for one-off sensitive queries you don't want kept in your regular history.
- Cowork gives Claude access to a specific folder on your computer. Point it at a scoped working directory rather than your whole file system, and check what it touched when a task finishes. By default, Cowork will request access to your user folder containing your Documents, Desktop, Downloads, and similar directories. Do not allow this.
- If you connect Cowork to Google Drive, Slack, or similar, remember that this provides live access to the data in those accounts going forward. You should instead provide limited access to where your public data are stored.

## Using Claude Code

Claude Code is the agentic coding tool, available as a terminal command-line interface (CLI), inside the Claude Desktop app's Code tab, or as an extension in Visual Studio Code and JetBrains.

### How to Install

- Follow the instructions on [Anthropic's website](https://claude.com/product/claude-code) for your operating system of interest.
- After installation, you can confirm the software has installed with `claude --version`.
- Your project directory should only include files you wish to send to Claude.
- You can run `claude` inside a project directory to start your first session.
- Claude Code requires a paid plan (not the free tier).
- If you'd rather not install the CLI separately, the Desktop app's Code tab gives you similar functionality with no extra setup.

### Protecting Your Data

Protecting your data in Claude Code is more technical.

Add explicit deny rules in `settings.json` so Claude Code can't read sensitive paths, regardless of what you or Claude ask for. To do this, open `settings.json` in the hidden `~/.claude/` directory and add the following JSON:

```json
{
  "cleanupPeriodDays": 7,
  "enableAllProjectMcpServers": false,
  "env": {
    "CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC": "1"
  },
  "permissions": {
    "defaultMode": "default",
    "deny": [
      "Read(./.env)",
      "Read(./.env.*)",
      "Read(./secrets/**)",
      "Read(~/.ssh/**)",
      "Read(~/.aws/**)"
    ],
    "disableBypassPermissionsMode": "disable"
  }
}
```

- You should also use environment variables to control how Claude Code behaves. You can disable non-essential outbound traffic with `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1`.
- `/feedback` sends your conversation history, including code, to Anthropic. Known API keys and tokens get redacted first, but everything else goes up as-is. Turn it off with `DISABLE_FEEDBACK_COMMAND=1` if you're working on anything you wouldn't want leaving your machine.
- Session transcripts get cached locally in plaintext under `~/.claude/projects/` for 30 days by default. If that's a concern for a given project, shorten `cleanupPeriodDays` or set `CLAUDE_CODE_SKIP_PROMPT_HISTORY`.
- Telemetry and error reporting only send operational metrics, not code or file paths, but you can turn those off too with `DISABLE_TELEMETRY=1` and `DISABLE_ERROR_REPORTING=1`.
- On Pro, Max, and Team plans, Claude Code starts sessions in auto mode, which uses a second model called the classifier to approve or block actions without prompting you. For research data we recommend switching to manual mode (press <kbd>Shift</kbd>+<kbd>Tab</kbd> to cycle modes, or set `permissions.defaultMode` in `settings.json`) to review each action yourself.

### Protecting Your API Key

If you authenticate Claude Code with an API key instead of a Claude.ai login, protect that key like a password.

- Never commit keys to a repository. Keep them out of source control. Use a config file like `settings.json`, environment variables, or a secrets manager.
- Never share your key or paste it into tickets, shared documents, or chat.
- Rotate your key periodically, and immediately if you suspect it has been exposed.
- Use one key per person or service where possible, so usage is attributable and a single key can be revoked without disrupting everyone else.
- Scrub keys from logs and screenshots before sharing them for troubleshooting. This matters especially given the session transcript caching and `/feedback` behavior described above.

[data-finder]: https://access.tufts.edu/data-finder
[genai-tufts]: https://it.tufts.edu/ai/generative-ai-tools
