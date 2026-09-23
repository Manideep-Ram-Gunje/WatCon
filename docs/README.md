# Documentation

Twenty-odd files. This says which one answers your question.

---

## I want to *use* it

| Read | For |
|---|---|
| [`../README.md`](../README.md) | what the tool does and what it found — start here |
| [`installation.rst`](installation.rst) | installing, including **PyMOL** and a clean-machine checklist |
| [`getting_started.rst`](getting_started.rst) | ten seconds to a result, then your own data |
| [`consurf_data.rst`](consurf_data.rst) | getting ConSurf results — the one manual step |
| [`faq/conservation_options.rst`](faq/conservation_options.rst) | every setting the extension adds, and what each refuses to guess |
| [`MANUAL_TESTING.md`](MANUAL_TESTING.md) | checking an installation really works, step by step |
| [`tutorials.rst`](tutorials.rst) | five worked tutorials |

**Shortest path:** install, run `watcon demo`, read
[`getting_started.rst`](getting_started.rst). About fifteen minutes.

---

## I want to *understand* it

| Read | For |
|---|---|
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | how the pieces fit, and the four rules behind the design |
| [`CODE_MAP.md`](CODE_MAP.md) | which file owns which capability — **the map to open first** |
| [`CONSURF_INTEGRATION.md`](CONSURF_INTEGRATION.md) | why each decision was made, at length |
| [`api.rst`](api.rst) | generated API reference for every module |
| [`../experiments/barnase_waters/FINDINGS.md`](../experiments/barnase_waters/FINDINGS.md) | the headline result, its method and its limits |
| [`../experiments/benchmark/FINDINGS.md`](../experiments/benchmark/FINDINGS.md) | PTP1B at scale, and the family of fifteen |

**Shortest path:** ARCHITECTURE, then CODE_MAP, then whichever module they send
you to. Half a day for the whole codebase.

---

## I want to *change* it

| Read | For |
|---|---|
| [`developer_guide.rst`](developer_guide.rst) | setup, tests, conventions, how to add things, what to be careful about |
| [`CODE_MAP.md`](CODE_MAP.md) | the "if you want to change X, touch Y" table |
| [`CONSURF_CHANGELOG.md`](CONSURF_CHANGELOG.md) | every change, why, and which alternative was rejected |

**Before your first change:** read the four invariants in ARCHITECTURE and the
warnings in the developer guide. Most of the defects in the changelog are one of
those rules being broken, and they fail *silently* — a wrong answer rather than
an error.

---

## Things worth knowing wherever you start

**This is a derivative work.** WatCon — the water-network analysis, clustering,
MSA machinery and PyMOL projections — is by Brownless, Harrison-Rawn and
Kamerlin (*JACS Au* 2025). The ConSurf integration is what this fork adds.
[`CODE_MAP.md`](CODE_MAP.md) marks which files are which. Cite the original.

**The findings include a negative result and a retraction**, both deliberately
kept: conservation does *not* help predict where water sits (ΔAUC ≈ 0), and a
claim that held at five proteins did not survive at fifteen. Neither is buried.

**Burial is not disentangled** in any result here, so no causal claim is made
anywhere. If you extend this work, that is the experiment worth doing.
