# Project Instructions

This repo contains portable AI-agent skills for clinical trial design R packages.

## Structure

The model-neutral source of truth is `.agents/skills/<skill-name>/` with:
- `SKILL.md` - Main skill definition (triggers, workflow, key concepts)
- `agents/openai.yaml` - Codex-compatible display metadata
- `references/` - Detailed code patterns, templates, and examples

Model-specific directories such as `.claude/skills/` are compatibility/export
targets, not the preferred place for new source edits.

VS Code-style agents should read skills directly from `.agents/skills/`.
Cursor should use `.cursor/rules/gsdesign-skills.mdc` as an adapter to those
same skills. Claude Code can use mirrored copies in `.claude/skills/`. Codex
can use the same `.agents/skills/` source and may also read
`agents/openai.yaml`.

## Cross-package skills

Skills that integrate multiple packages use hyphenated names (e.g., `graphicalMCP-gsDesign2`).

## Adding a new skill

1. Create a directory under `.agents/skills/`
2. Add a `SKILL.md` with YAML frontmatter (`name`, `description`)
3. Add `agents/openai.yaml` for display metadata
4. Add `references/` subdirectory for code patterns when needed
5. Update `README.md`
