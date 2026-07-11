# gsDesignSkills

AI skills for clinical trial design R packages:

- **gsDesign** - Group sequential design (classical)
- **gsDesign2** - Group sequential design (next generation)
- **rpact** - Confirmatory adaptive clinical trial design and analysis
- **graphicalMCP** - Graphical multiple comparison procedures
- **simtrial** - Clinical trial simulation
- **wpgsd** - Weighted parametric group sequential design
- **gMCPLite** - Graphical MCP (legacy)
- **illness-death** - Illness-death multi-state model for correlated endpoints
- **gsDesignNB** - Negative binomial / recurrent event designs

## Cross-package skills

- **graphicalMCP-gsDesign2** - Group sequential design with graphical multiplicity control (Maurer-Bretz framework)
- **multi-endpoint-sim** - Multi-endpoint trial simulation with sequential testing (gsDesign + illness-death + simtrial + graphicalMCP)

## Semantic layer

- **semantic-router** - Maps natural-language trial design requests to normalized intent, package skills, and workflows
- **trial-design-glossary** - Defines shared endpoint, estimand, timing, information, multiplicity, adaptation, and simulation concepts

The semantic layer also includes:

- `glossary/` - Shared trial-design vocabulary
- `crosswalks/` - Package-specific concept-to-function mappings
- `evals/` - Natural-language routing examples and expected outputs

## Usage

Use `.agents/skills/` as the model-neutral source of truth. Each skill keeps
portable Markdown instructions in `SKILL.md` and optional references, plus
`agents/openai.yaml` metadata for Codex UI surfaces.

Different assistants can consume the same source in slightly different ways:

- **VS Code agents**: read `AGENTS.md` and the relevant `.agents/skills/<skill-name>/SKILL.md`.
- **Cursor**: use `.cursor/rules/gsdesign-skills.mdc`, which routes Cursor to `.agents/skills/`.
- **Claude Code**: copy or mirror selected skills into `.claude/skills/`.
- **Codex**: use `.agents/skills/` directly, including `agents/openai.yaml` for display metadata.

```bash
# Copy a single skill for a Claude Code project
cp -r .agents/skills/graphicalMCP-gsDesign2 /path/to/your/project/.claude/skills/

# Copy all portable skills
cp -r .agents/skills/* /path/to/your/project/.claude/skills/
```

Then reference the skill in your project's `CLAUDE.md`:

```markdown
## Skills

When working on group sequential designs with graphical multiplicity control,
read the skill at `.claude/skills/graphicalMCP-gsDesign2/SKILL.md`.
```

## API documentation

Several skills include vendored `references/llms.txt` API references, plus
local `llms_local.txt` files where package docs have been regenerated. Some
older `llms.txt` files originated from gsDesign.ai, but gsDesignSkills is the
source of truth for the skills and does not require that external site.

## Contributing

Add new portable skills by creating a directory under `.agents/skills/` with a `SKILL.md` file.
Use a `references/` subdirectory for detailed code patterns and examples.
Add `agents/openai.yaml` with display metadata for Codex-compatible surfaces.
If a model-specific copy is needed, export from `.agents/skills/` after review.
Update `.cursor/rules/gsdesign-skills.mdc` only when Cursor-specific routing or loading behavior changes.
