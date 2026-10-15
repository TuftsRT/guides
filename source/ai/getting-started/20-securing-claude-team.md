---
tags: ai generative-ai claude claude-code data-privacy administration
---

# Securing Claude Team for Use at Tufts

This document is for Tufts faculty bringing graduate students, postdocs, or research staff onto Claude Enterprise (Claude Team). Claude Enterprise refers to the paid Claude Team plan.

This does not refer to the Pro and Max plans, which are consumer-level and do not ensure your data are not used for training or product development. We recommend using Claude Team for research groups at Tufts.

For general guidance on similar tools, see Using Generative AI and Agentic AI in Research. This document supplements the university-wide [Tufts Generative AI Usage Guidelines for Faculty, Staff, and Students](https://ovpe.tufts.edu/initiatives/ai/ai-guidelines/) and [Generative AI Guidance](https://it.tufts.edu/ai) with recommended settings from Tufts Technology Services (TTS).

```{note}
This document assumes Owner or Admin access to your organization's settings. If you don't have it, or you have any questions, contact Research Technology at <tts-research@tufts.edu>.
```

## Getting Started

Five key topics are critical when you are setting up your Claude Team for the first time. Each has its own section below.

- **Sharing:** Restrict sharing, so your team's project or chat can't end up visible outside your lab.
- **Spending:** Set spend limits by group, so a runaway agentic task doesn't run up an open-ended bill.
- **Training:** Confirm training is off and set a retention period, so you know what's stored, for how long, and who can see it.
- **Coding:** Lock down Claude Code before first sign-in, with a policy the users can't override.
- **Connecting:** Default connectors to read-only and know what data each one reaches. Allow users to request connectors, but do not blanket allow.

## Restrict Sharing

- Turn off [Public projects](https://support.claude.com/en/articles/9927533-disable-public-projects-for-your-organization) under **Organization settings > Data and Privacy**. This converts existing public projects to private and stops anyone from creating new ones. Members can still share a private project with specific people, so this doesn't block normal lab collaboration.
- If you create a project for a student's research, set its [visibility to Private](https://support.claude.com/en/articles/9519189-manage-project-visibility-and-sharing) and invite the student directly. Before anything gets uploaded, check the dataset's classification against the [Data Storage Finder][data-finder] and follow all regulations.
- Turn off [Rate chats](https://support.claude.com/en/articles/10504844-managing-user-feedback-settings-on-claude-for-work-team-and-enterprise-plans) in the same settings if you don't want students sending full chat content to Anthropic through the thumbs up / thumbs down feedback button.

## Set Spend Limits by Group

Spend limits cascade: organization, then group, then individual. The lowest limit wins, so an organization-wide cap protects you even if a group or user limit gets set higher by mistake.

- Create a group for the cohort under **Organization settings > Groups** and set a per-user monthly limit for it under **Organization settings > Usage > By group**. See [Manage groups and group spend limits on Enterprise plans](https://support.claude.com/en/articles/13799932-manage-groups-and-group-spend-limits-on-enterprise-plans). This covers every student you add going forward, so you're not setting limits one at a time.
- Start conservative. It's easier to raise a student's limit on request than to walk back an unexpected charge. Override the group limit for a specific student under the **By member** tab if their project needs more.
- If students will use Claude Code heavily, its per-user API key lives in its own [workspace](https://platform.claude.com/docs/en/manage-claude/workspaces), and that workspace's share of the organization's limits can be capped separately under **Settings > Workspaces**.

## Training and Retention

Claude for Work plans, including Enterprise, don't train on your prompts, files, or Claude's replies by default. See [Is my data used for model training?](https://support.anthropic.com/en/articles/7996868) Confirm nothing's been toggled at the organization level. The exception is feedback: a thumbs up / thumbs down or bug report sends the full conversation to Anthropic, which is what Rate chats (above) turns off.

- Set a [custom retention period](https://privacy.anthropic.com/en/articles/10440198-custom-data-retention-controls-for-claude-enterprise) under **Organization settings > Data and Privacy**. Training and retention are separate settings, and the default retention period isn't automatically right for a research cohort.
- Plan for offboarding before you need it. Individual members can't self-serve export their data; only your organization's Primary Owner can [export organization data](https://support.claude.com/en/articles/13346720-export-your-organization-s-data). Add a step to your lab's offboarding checklist to pull a graduating student's chats and projects, and to check for anything they made public, before you remove their account.

Tell students plainly what Tufts' role as Primary Owner means: the same kind of access the university holds over Tufts email or Drive. Tufts can retrieve data for public safety or litigation reasons, but doesn't routinely review conversation content.

## Locking Down Claude Code

The companion guide shows a student how to edit their own `settings.json` under `~/.claude/`, which they can change or remove at any time. For a cohort, deploy the same restrictions as an [enterprise managed policy](https://code.claude.com/docs/en/settings) instead, so they're in place before the student's first session and can't be overridden.

Deploy `managed-settings.json` to the system path for your platform, using whatever device management Tufts already uses (mobile device management, Ansible, or Group Policy Objects):

- macOS: `/Library/Application Support/ClaudeCode/managed-settings.json`
- Linux and Windows Subsystem for Linux (WSL): `/etc/claude-code/managed-settings.json`
- Windows: `C:\Program Files\ClaudeCode\managed-settings.json`

```{warning}
On Windows, use `C:\Program Files\ClaudeCode\`. Claude Code does not read the legacy `C:\ProgramData\ClaudeCode\managed-settings.json` path. A policy deployed there fails silently: students get no restrictions, and nothing reports an error.
```

A starting policy for a research environment:

```json
{
  "cleanupPeriodDays": 7,
  "allowManagedMcpServersOnly": true,
  "enableAllProjectMcpServers": false,
  "allowManagedPermissionRulesOnly": true,
  "env": {
    "CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC": "1",
    "DISABLE_TELEMETRY": "1"
  },
  "permissions": {
    "defaultMode": "default",
    "deny": [
      "Read(./.env)",
      "Read(./.env.*)",
      "Read(./secrets/**)",
      "Read(~/.ssh/**)",
      "Read(~/.aws/**)",
      "Bash(cat .env*)",
      "Bash(printenv)",
      "Bash(env)",
      "Bash(grep -r AKIA*)"
    ],
    "disableBypassPermissionsMode": "disable",
    "disableAutoMode": "disable"
  }
}
```

- `allowManagedMcpServersOnly` and `allowManagedPermissionRulesOnly` stop a student from adding their own Model Context Protocol (MCP) servers or permission rules that would override the deny list.
- `disableBypassPermissionsMode: "disable"` prevents anyone from starting a session with `--dangerously-skip-permissions`, on this machine or any project.
- `disableAutoMode: "disable"` removes auto mode entirely, so no one in the cohort can select it. On Pro, Max, and Team plans auto mode is the starting mode, and a classifier approves actions instead of the student. The companion guide tells students to switch to manual mode themselves, but that choice is theirs to undo; this key settles it centrally. Use `permissions.defaultMode` instead if you want students to start in manual mode while still being able to switch.
- A `Read` deny rule only blocks Claude Code's Read tool. It doesn't stop the same file coming back through Bash, so the deny list also covers `cat`, `printenv`, `env`, and `grep` against the same patterns. Test any deny rule with dummy credentials before trusting it with real ones.
- Add paths specific to your lab, such as a directory holding restricted datasets, and keep `.env` files and other secrets outside the directory Claude Code is scoped to where you can. A secrets manager is safer still: nothing sensitive sits on disk to be read in the first place.

```{important}
Verify the policy took effect. Run `claude` and check `/status` inside a session; it should show the managed settings source before you tell a student to start working. Remember, Claude Code can reach outside the working directory folder easily. Be sure to read every line of code it creates, as it is easy to mistakenly expose private datasets.
```

## Connectors

Connectors (Google Drive, Slack, and similar; see [Use connectors to extend Claude's capabilities](https://support.claude.com/en/articles/11176164-use-connectors-to-extend-claude-s-capabilities)) inherit whatever access the connected account already has, and aren't covered by anything above. Default them to read-only. A connector that can send messages or write files as the student is an exfiltration path if a session gets manipulated by a hidden instruction in a document or message it reads.

- Start narrow. Before connecting a shared Drive or Slack workspace, know what's actually in it, and connect a specific folder rather than root-level access where the connector allows it.

## Other Things to Know

Items that come up less often:

- Research and citation connectors (PubMed, bioRxiv, and similar) carry a distinct risk: the query text itself, not just retrieved documents, is sent to and may be logged by that vendor. Unpublished results or restricted data don't belong in the search box either. Use the [local AI models](https://go.tufts.edu/AITools) for research you don't want to release yet.
- A connector's general risk approval doesn't automatically cover protected health information (PHI). If a student's project touches health data, confirm that separately.
- A spend limit caps dollars, but Claude also enforces its own usage limits underneath (hourly, daily, weekly, monthly). If a student hits one, the only fix is to wait or pay more for extra credits in advance.
- If a cohort's data requires Zero Data Retention or HIPAA-equivalent terms, publishing a chat as a shared Artifact link won't be available for that group.
- Don't assume every Claude product or model is available just because you're on Enterprise.

## About Chat (Claude.ai Web Portal)

A few things live in the chat interface itself, not the admin console:

- Use [Incognito chat](https://support.claude.com/en/articles/12260368-use-incognito-chats) (the ghost icon, or <kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>I</kbd>) for one-off questions involving anything sensitive. It's never used to improve Claude, but Anthropic still retains it for a safety-review window, so it supplements good judgment rather than replacing it.
- If memory is enabled, turn off ["Include sensitive topics in memory"](https://support.claude.com/en/articles/11817273-use-claude-s-chat-search-and-memory-to-build-on-previous-context) under **Settings > Memory**, and check [**Settings > Privacy > Manage Shared chats**](https://support.claude.com/en/articles/10593882-share-and-unshare-chats) occasionally for anything you don't want reachable by a link.
- Redact or replace names, identifiers, and account numbers with placeholders before pasting a document into a chat.

## Delegating Admin Access

Give a co-principal investigator (co-PI) or lab manager only the access their role needs. [Custom roles](https://support.claude.com/en/articles/13930452-manage-custom-roles-on-enterprise-plans) grant all View or all Manage within an area (Privacy, User Management, Libraries, and so on), so build the role deliberately rather than defaulting to broad Admin access.

## Questions

For help with any of the settings above, or with questions before using Claude with research data, contact <tts-research@tufts.edu>.

[data-finder]: https://access.tufts.edu/data-finder
