# Using the Trial Design Glossary

Read `../../../glossary/core-concepts.md` when a request uses domain terms that
could map to more than one package or workflow.

## Normalization principles

- Prefer statistical intent over package-specific argument names.
- Keep the user's clinical terms when they carry meaning, then add the
  normalized semantic field.
- Treat synonyms as routing cues, not as exact substitutes. For example,
  "events" can route to survival or recurrent-event workflows depending on the
  rest of the request.
- Ask a clarification question when ambiguity changes the endpoint family,
  estimand, multiplicity method, adaptation method, or timing basis.

## Common synonym cues

- Recurrent events, exacerbations, hospitalizations, counts over follow-up:
  `recurrent_event` or `negative_binomial_count`.
- Time to event, survival, OS, PFS, hazard ratio: `time_to_event`.
- Graph, nodes, transition weights, alpha recycling: `graphical_mcp`.
- Blinded information, pooled rate, nuisance re-estimation:
  `information_source: blinded`.
- Information fraction, spending time, interim look: group sequential timing.
