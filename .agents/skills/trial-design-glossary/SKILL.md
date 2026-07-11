---
name: trial-design-glossary
description: >
  Define and normalize shared clinical trial design vocabulary for gsDesign
  ecosystem skills. Use when Codex needs consistent meanings for endpoint
  family, estimand, hypothesis type, information, timing, boundaries,
  multiplicity, adaptation, simulation, or other semantic fields before
  choosing a package workflow.
---

# Trial Design Glossary

Use this skill to standardize terminology before routing or implementation.

## Resources

- Canonical glossary: `../../../glossary/core-concepts.md`
- Usage notes: `references/using-glossary.md`

## Workflow

1. Map user terms to the closest canonical concept.
2. Preserve clinically meaningful details such as endpoint family, estimand,
   time scale, population, and intercurrent-event handling.
3. Flag terms that are ambiguous enough to change the package or workflow.
4. Hand off to `semantic-router` when package selection is needed.

Do not perform package-specific calculations from this skill. Use the package
skill selected by the router.
