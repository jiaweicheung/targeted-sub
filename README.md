# Buyer-Specific Payments Separate Equilibrium Convexity from Posted-Price Integrability

Working paper, 38 pages. This folder holds the full LaTeX source, the
bibliography, the compiled PDF, and the exact-arithmetic verification scripts.

## Build

Main file: `convex_potentials_to_Jiawei.tex`

```text
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
bibtex convex_potentials_to_Jiawei
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
```

The generated `.bbl` is included, so the build also works without running
BibTeX. Substantive text lives entirely in `convex_potentials_body.tex`; the
wrapper and `preamble_common.tex` only carry the class options and macros.

## Files

| File | Role |
|---|---|
| `convex_potentials_to_Jiawei.tex` | wrapper, inputs the preamble and the body |
| `preamble_common.tex` | document class, packages, macros, theorem environments |
| `convex_potentials_body.tex` | the whole manuscript |
| `convex_potentials.bib` | 32 references, each with a DOI |
| `convex_potentials_to_Jiawei.bbl` | generated bibliography |
| `convex_potentials_to_Jiawei.pdf` | compiled manuscript |
| `README_verification.md` | what each verification script checks |

## Verification scripts

Five scripts check the finite examples and the appendix robustness lemma in
exact rational arithmetic. They support the proofs and do not replace them; the
convex-analytic steps remain analytic. See `README_verification.md` for what
each one covers and how to run them.
