# gsDesignSkills

Claude Code skills for clinical trial design R packages:

- **gsDesign** - Group sequential design (classical)
- **gsDesign2** - Group sequential design (next generation)
- **rpact** - Confirmatory adaptive clinical trial design and analysis
- **graphicalMCP** - Graphical multiple comparison procedures
- **simtrial** - Clinical trial simulation
- **wpgsd** - Weighted parametric group sequential design
- **gMCPLite** - Graphical MCP (legacy)

## Cross-package skills

- **graphicalMCP-gsDesign2** - Group sequential design with graphical multiplicity control (Maurer-Bretz framework)

## Usage

Copy the skills you need into your project's `.claude/skills/` directory:

```bash
# Copy a single skill
cp -r .claude/skills/graphicalMCP-gsDesign2 /path/to/your/project/.claude/skills/

# Copy all skills
cp -r .claude/skills/* /path/to/your/project/.claude/skills/
```

Then reference the skill in your project's `CLAUDE.md`:

```markdown
## Skills

When working on group sequential designs with graphical multiplicity control,
read the skill at `.claude/skills/graphicalMCP-gsDesign2/SKILL.md`.
```

## Contributing

Add new skills by creating a directory under `.claude/skills/` with a `SKILL.md` file.
Use a `references/` subdirectory for detailed code patterns and examples.
