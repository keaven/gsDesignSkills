---
name: semantic-router
description: >
  Route natural-language clinical trial design requests into normalized
  semantic intent, package-specific skills, and recommended gsDesign ecosystem
  workflows. Use when a request spans multiple trial design concepts, when the
  right package is unclear, or when Codex needs to map endpoint family,
  estimand, information, timing, boundaries, multiplicity, adaptation, or
  simulation requirements to gsDesign, gsDesign2, gsDesignNB, simtrial,
  graphicalMCP, rpact, wpgsd, illness-death, or related skills.
---

# Trial Design Semantic Router

Use this skill before package-specific skills when user intent needs to be
translated into a trial-design workflow.

## Resources

- Routing rules: `references/routing-patterns.md`
- Shared vocabulary: `../../../glossary/core-concepts.md`
- Package crosswalks: `../../../crosswalks/*.yaml`
- Routing examples: `../../../evals/natural-language-requests.yaml`

## Workflow

1. Normalize the request into semantic fields: endpoint family, estimand,
   hypothesis type, design family, timing basis, information source,
   boundaries, multiplicity, adaptation, and simulation.
2. Read the relevant glossary entries and package crosswalks only as needed.
3. Select the smallest set of package skills that can satisfy the intent.
4. Identify missing inputs that change the workflow, analysis, or package
   choice. Ask only for those inputs before implementation.
5. Hand off to the selected package skill for code patterns and API details.
6. Explain the routing decision in terms of user intent and package capability.

## Output Pattern

When useful, summarize routing with this compact structure:

```yaml
intent:
  endpoint_family:
  estimand:
  design_family:
  timing_basis:
  information_source:
  multiplicity:
  adaptation:
  simulation:
routing:
  skills:
  packages:
  workflow:
  missing_inputs:
```

Keep the structure short. Do not force every field when the user asked a narrow
question.

## Clarification Rule

Ask a clarification question only when the answer changes the recommended
package, workflow, statistical method, or required inputs. Otherwise, make a
reasonable assumption and state it.
