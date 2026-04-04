# Project Instructions

This repo contains Claude Code skills for clinical trial design R packages.

## Structure

Each skill lives in `.claude/skills/<skill-name>/` with:
- `SKILL.md` - Main skill definition (triggers, workflow, key concepts)
- `references/` - Detailed code patterns, templates, and examples

## Cross-package skills

Skills that integrate multiple packages use hyphenated names (e.g., `graphicalMCP-gsDesign2`).

## Adding a new skill

1. Create a directory under `.claude/skills/`
2. Add a `SKILL.md` with YAML frontmatter (`name`, `description`)
3. Add `references/` subdirectory for code patterns
4. Update `README.md`
