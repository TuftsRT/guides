---
tags: ai api openai claudecode assistant hpc agent security
---

# AI Coding Assistants on Tufts HPC

**Author:** Kyle Monahan, Research Technology, TTS

AI coding assistants (sometimes called "agentic coding tools" or "AI pair programmers") can meaningfully accelerate research workflows on the Tufts HPC cluster, helping you write Slurm scripts, debug Python pipelines, structure batch jobs, and navigate the Linux environment. This guide covers when and how to use them responsibly on the Tufts HPC.

```{warning} AI Usage Disclaimer
AI-generated content may be incorrect, biased, incomplete, or misleading. Always review, verify, and edit AI output before use. Before uploading data to any AI tool, please review the [Tufts Information Classification and Handling Policy](https://it.tufts.edu/information-classification-and-handling) and the [Data Storage Finder](https://it.tufts.edu/data-storage-finder). Do not enter any restricted data into the tool without formal TTS review. For requirements and best practices, please refer to Tufts' [Guidelines for Use of Generative AI Tools](https://it.tufts.edu/guidelines-use-generative-ai-tools).
```

---

## Do I Need an AI Coding Assistant?

Not everyone needs an AI Coding Assistant. AI coding tools are most valuable when you:

- Are new to HPC, Linux, or scripting and want guided help writing and understanding code
- Need to rapidly prototype or adapt an existing script to run on the cluster
- Are writing repetitive boilerplate (Slurm headers, argument parsing, file I/O loops)
- Want to debug error messages or understand unfamiliar software output
- Are translating a workflow from one language or framework to another

If you are an experienced HPC user with established scripts and workflows, you may find these tools add less value, and can even slow you down, as you will need to spend time reviewing and correcting AI-generated code.

**The key principle:** AI assistants are drafting tools, not authoritative sources. **Always review generated code before running it on the cluster, especially for anything involving file paths, data handling, or external API calls.** They are not a replacement for your critical thought and rigor.

---

## Recommended Tools

There is no single "best" tool. The right choice depends on your workflow, editor preferences, and budget. The following are the most widely used options in research settings as of 2025–2026.

### GitHub Copilot *(Most accessible; recommended starting point)*

GitHub Copilot integrates directly into popular editors including VS Code, JetBrains IDEs (PyCharm, IntelliJ), Neovim, and Visual Studio. It provides inline code suggestions and a conversational chat panel without requiring you to change editors or workflows.

- **Best for:** Researchers already using VS Code or JetBrains who want low-friction AI assistance; teams with mixed editor preferences
- **Pricing:** Free tier available to all GitHub users (2,000 completions/month, 50 chat messages/month); Pro plan ~\$10/month; **free for verified students and educators via [GitHub Education](https://education.github.com/)**.
- **Data note:** By default, Copilot Free/Pro may use interaction data for model improvement. Check your [GitHub settings](https://github.com/settings/copilot) and opt out of telemetry if working with sensitive data. Copilot Business/Enterprise plans have stronger data protections.
- **Getting started:** [GitHub Copilot Docs](https://docs.github.com/en/copilot)
- **Microsoft Copilot note:** GitHub Copilot is a coding assistant and is distinct from Microsoft Copilot (M365). All Tufts community members have access to basic Microsoft Copilot Chat through their existing Microsoft license - see the [Tufts Microsoft Copilot guide](https://it.tufts.edu/guides/microsoft-copilot).
- **HPC support:** Not formally supported on the Tufts HPC cluster, though you can configure it on your own workstation.

### Cursor *(AI-native IDE)*

Cursor is a standalone editor built as a VS Code fork, redesigned around AI assistance. Its "Composer" feature can plan and edit across multiple files simultaneously - useful when refactoring or adapting complex multi-script workflows.

- **Best for:** Researchers with larger, multi-file codebases who want deep AI integration and don't mind adopting a new editor
- **Pricing:** Free tier; Pro ~\$20/month
- **Data note:** Cursor is SOC 2 Type II compliant. Review their [privacy policy](https://cursor.com/privacy) before use with sensitive research data.
- **Getting started:** [Cursor Docs](https://docs.cursor.com)
- **HPC support:** Not formally supported on the Tufts HPC cluster, though you can configure it on your own workstation.

### Windsurf *(AI-native IDE; best onboarding for newcomers)*

Windsurf is an AI-native IDE built from scratch by Codeium around its **Cascade** feature, which enables multi-step AI workflows with real-time feedback. Unlike Cursor (a VS Code fork), Windsurf was built independently, and is widely noted for having the most polished onboarding experience among AI coding tools - making it a strong choice for researchers who are new to AI-assisted development.

- **Best for:** Researchers new to AI coding tools; rapid prototyping; projects where you want the AI to watch for errors and iterate automatically without manual copy-pasting
- **Pricing:** Free tier with generous limits; Pro ~\$15/month
- **Data note:** Review [Codeium's privacy policy](https://codeium.com/privacy-policy) before use with sensitive research data.
- **Getting started:** [Windsurf Docs](https://docs.codeium.com/windsurf/getting-started)
- **HPC support:** Not formally supported on the Tufts HPC cluster, though you can configure it on your own workstation.

### Claude (claude.ai) and ChatGPT *(Best for conversational help and script explanation)*

Web-based chat interfaces like [claude.ai](https://claude.ai) and [ChatGPT](https://chatgpt.com) are excellent for asking questions in plain English - explaining error messages, generating starter scripts, or discussing HPC concepts. They require no installation and work from any browser.

- **Best for:** Getting quick explanations, learning HPC concepts, drafting initial scripts before testing on the cluster, troubleshooting errors by pasting output
- **Pricing:** Free tiers available; paid plans unlock more capable models
- **Data note:** Do **not** paste sensitive data, HIPAA-covered data, unpublished research data, or proprietary code into public-facing chat interfaces. See [Data and Privacy](#data-and-privacy) below.

### Claude Code *(Advanced; terminal-native agentic coding)*

Claude Code is a command-line tool from Anthropic that operates directly in your terminal, can read your entire codebase, and autonomously edits files across a project. It is the most capable option for complex, multi-file tasks, but requires familiarity with terminal workflows and has usage-based pricing.

- **Best for:** Advanced users comfortable in the terminal who need autonomous multi-file editing or codebase-wide refactoring
- **Pricing:** Usage-based via Anthropic API (~$3–$15/session depending on codebase size); Claude Pro subscription (~\$20/month) provides included usage
- **Note:** Claude Code runs locally on your workstation or laptop - **not** on the HPC cluster itself. Use it to write and refine scripts before transferring them to the cluster.
- **Getting started:** [Claude Code Docs](https://docs.anthropic.com/en/docs/claude-code)

---

## How to Use AI Assistants with the HPC Cluster

### Recommended workflow

1. **Draft locally with AI assistance.** Use your AI tool of choice to write or adapt your script on your laptop.
1. **Review the generated code carefully.** Check file paths, resource requests, and any logic involving your data. AI tools can produce plausible-looking but incorrect code.
1. **Test with a small job first.** Transfer the script and run it on a single file or small subset before scaling up to a full array.
1. **Iterate.** If the job fails, paste the error output back into your AI assistant for help debugging.
1. **Scale up.** Once the single job succeeds, submit your full job array.

### Useful prompts for HPC work

AI assistants respond well to specific, contextual prompts. Some examples:

- *"Write a Slurm batch script for a Python job that processes one file per task in a job array. Each task needs 1 CPU, 4GB memory, and up to 2 hours."*
- *"My Python script exits with `OSError: [Errno 28] No space left on device`. What does this mean on an HPC cluster and how do I fix it?"*
- *"Convert this for-loop that processes 500 files sequentially into a Slurm job array."*
- *"Explain what `#SBATCH --array=0-499%20` does."*
- *"I need to store my API key securely in a Slurm job. What is the best practice?"*

---

## Data and Privacy

### What you should NOT share with public AI tools

- **HIPAA-covered data** (patient records, clinical data, identifiable health information) - never input into any public AI tool
- **IRB-protected human subjects data** - check your IRB protocol; sharing with third-party AI may violate approval conditions
- **Unpublished research data or manuscripts** - may affect intellectual property rights or publication priority
- **Proprietary datasets** covered by data use agreements or sponsor restrictions
- **Passwords, API keys, or credentials** of any kind

### What is generally safe to share

- Generic code logic, Slurm script structure, and programming questions
- Anonymized or synthetic data used for testing
- Publicly available datasets you are already permitted to use
- Error messages and log output (check these don't contain file paths revealing sensitive information)

### A practical rule of thumb

> **If you would not paste it into a public GitHub repository, do not paste it into a public AI chat interface.**

If you work with sensitive or regulated data, consult [Tufts Research Technology](mailto:tts-research@tufts.edu) about appropriate options, including enterprise or self-hosted AI tools with stronger data protections.

---

## Common Pitfalls

**AI-generated Slurm scripts may request wrong resources.** AI tools are trained on general examples and may suggest excessive memory, CPU counts, or wall times. Requesting far more than you need wastes your fairshare allocation and delays queue start times. Always benchmark with small runs first.

**AI-generated paths may not exist on the cluster.** Tools like Copilot and ChatGPT will sometimes invent plausible-looking paths (`/scratch/user/`, `/home/user/data/`) that do not match Tufts HPC's actual directory structure. Verify paths against the Storage documentation.

**AI assistants can hallucinate module names.** When asking for help loading software with `module load`, AI tools may suggest module names that do not exist on Tufts HPC. Use `module spider <software>` on the cluster to verify availability.

**Generated code may introduce subtle bugs.** Studies have found that AI-generated code can contain logic errors, insecure practices, or code duplication that is easy to miss during review. Treat AI output as a starting draft, not finished code.

**Do not run AI-suggested commands without understanding them.** If an AI suggests a command like `rm -rf`, `chmod 777`, or a complex pipeline you don't recognize, look it up before running it on the cluster.

---

## Summary: Quick Reference

| Tool                | Type       | Best For                          | Free Tier       | Key Consideration            |
| ------------------- | ---------- | --------------------------------- | --------------- | ---------------------------- |
| GitHub Copilot      | IDE plugin | VS Code / JetBrains users         | Yes (all users) | Review telemetry settings    |
| Cursor              | AI IDE     | Multi-file projects               | Yes             | Requires adopting new editor |
| Windsurf            | AI IDE     | Newcomers; rapid prototyping      | Yes (generous)  | Best onboarding experience   |
| claude.ai / ChatGPT | Web chat   | Q&A, script drafting, debugging   | Yes             | No sensitive data            |
| Claude Code         | CLI agent  | Terminal-first, complex codebases | No (API usage)  | Local only; advanced users   |

---

## Getting Help

If you have questions about using AI coding tools with the Tufts HPC cluster, or need guidance on whether a specific workflow or dataset is appropriate for use with AI tools:

- **Email:** [tts-research@tufts.edu](mailto:tts-research@tufts.edu)
- **Tufts AI Policy resources:** [Tufts Office of the Provost AI guidance](https://provost.tufts.edu)
- **Tufts information security:** [it.tufts.edu](https://it.tufts.edu)

> Linked external tools and resources are not affiliated with or endorsed by Tufts University. AI tools and their pricing, features, and data policies change frequently - verify current details with each vendor before use.
