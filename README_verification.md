# Computational verification

Five scripts perform exact-rational, symbolic, and seeded randomized checks for
the examples and algebra used in the manuscript. These checks support, but do
not replace, the proofs. The convex-analytic steps, namely facet-to-ridge
endpoint transport, convex gluing, and the closed-graph arguments, remain
analytic and are not delegated to code.

`verification.py` covers the smooth classification: the exact linear and kink
examples, the difference between a scaled price and the integrated effective
price at a kink, and randomized kink configurations and rank-one marginal
matrices. It writes `verification_report.txt`.

`verify_three_facet_certificate.py` covers the three-good log-affine
obstruction: all three ridge demand correspondences and excess-supply vertices,
the open ridge-to-ridge segments through affine endpoint margins, four
full-dimensional cell witnesses, normals derived from the stated $K$ matrix, all
six zero-coordinate ray isolations, the three multiplier ratios, and the cycle
products $32$ against $75$.

`verify_two_good_obstruction.py` covers the displayed three-buyer, two-good
market: the exchange axiom, exact demands and inclusion-minimal demands on every
path regime, the complete direction classification at the two reset prices, the
three crossing cells, and the one-dimensional excess-supply hull.

`verify_two_good_robustness.py` covers the appendix robustness lemma. It
recomputes the perturbed event prices from their closed forms, recovers the
least equilibrium price independently by vertex enumeration over the
assignment-supporting polyhedra, checks the map from update directions to demand
profiles by taking an exact small step, and cross-tabulates conditions (R1) to
(R5) against the conclusions over a perturbation battery.

`independent_math_audit.py` re-implements the finite checks separately as a
regression test against the other scripts.

## Requirements

- Python 3.10 or later
- NumPy and SymPy, for `verification.py` only

The other four scripts use the standard library alone, with exact arithmetic
through `fractions.Fraction`.

## Run

```text
python verification.py
python verify_three_facet_certificate.py
python verify_two_good_obstruction.py
python verify_two_good_robustness.py
python independent_math_audit.py
```

Randomized checks use fixed seeds, `20260817` in `verification.py` and
`20260902` in `verify_two_good_robustness.py`.
