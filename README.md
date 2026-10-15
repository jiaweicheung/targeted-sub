# Buyer-Specific Payments Separate Equilibrium Convexity from Posted-Price Integrability

This folder holds the full LaTeX source, the
bibliography, the compiled PDF, and the exact-arithmetic verification scripts.

## Build

Template-formatted working paper: `working-paper.tex`. It preserves the original
manuscript's content with the supplied template's 11-point A4 article layout,
Palatino text and mathematics, 1.2-inch margins, and one-and-a-half line spacing.

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error working-paper.tex
```

This produces `working-paper.pdf`. LaTeX intermediates are ignored by Git;
the `.tex` source and finished `.pdf` remain trackable. The existing tracked
files, including the original `.bbl` and PDF, are retained.

Original Elsevier-format manuscript: `convex_potentials_to_Jiawei.tex`

```text
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
bibtex convex_potentials_to_Jiawei
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
xelatex -interaction=nonstopmode -halt-on-error convex_potentials_to_Jiawei.tex
```

The original manuscript's generated `.bbl` is included, so that build also works
without running BibTeX. Its substantive text lives in
`convex_potentials_body.tex`; the wrapper and `preamble_common.tex` only carry
the class options and macros. `working-paper.tex` contains the same manuscript
in the new template and shares `convex_potentials.bib`. Keep substantive edits
in both manuscript sources in sync.

## Files

| File | Role |
|---|---|
| `working-paper.tex` | complete manuscript in the supplied working-paper template |
| `working-paper.pdf` | compiled template-formatted manuscript |
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
