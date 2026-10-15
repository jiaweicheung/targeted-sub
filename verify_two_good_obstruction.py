#!/usr/bin/env python3
"""Exact verification for the two-good lower-box and path obstruction.

This script independently checks the finite-market certificate used in R6:
  * model compatibility, including the M-natural exchange axiom;
  * exact demands, minimal demands, and inclusion-minimal maximally
    over-demanded sets on all path regimes;
  * the complete direction-regime classification at the two reset prices;
  * the three crossing cells and the one-dimensional excess-supply hull;
  * an optional exact consistency check for Algorithm 2's dual direction.

The minimum-equilibrium lower bound and the convex-analytic gluing/closed-graph
arguments are proved analytically in the manuscript.  They are not delegated to
this script.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product
from typing import FrozenSet, Iterable, Sequence, Tuple

Bundle = Tuple[int, int]
Price = Tuple[F, F]
Subset = FrozenSet[int]

BUNDLES: tuple[Bundle, ...] = ((0, 0), (1, 0), (0, 1), (1, 1))
SUPPLY: Bundle = (1, 1)
VALUES: tuple[tuple[int, int], ...] = ((10, 10), (10, 11), (0, 5))
SLOPES: tuple[tuple[int, int], ...] = ((1, 1), (1, 2), (1, 1))
SUBSETS: tuple[Subset, ...] = tuple(
    frozenset(j for j in range(2) if mask & (1 << j)) for mask in range(4)
)


def valuation(i: int, x: Bundle) -> F:
    if x == (0, 0):
        return F(0)
    return max(F(VALUES[i][j]) for j in range(2) if x[j])


def utility(i: int, x: Bundle, p: Price) -> F:
    return valuation(i, x) - sum(F(SLOPES[i][j]) * p[j] * x[j] for j in range(2))


def demand(i: int, p: Price) -> FrozenSet[Bundle]:
    values = {x: utility(i, x, p) for x in BUNDLES}
    best = max(values.values())
    return frozenset(x for x, u in values.items() if u == best)


def minimal_bundles(ds: Iterable[Bundle]) -> FrozenSet[Bundle]:
    items = tuple(ds)
    return frozenset(
        x for x in items
        if not any(y != x and all(y[j] <= x[j] for j in range(2)) for y in items)
    )


def mu(i: int, subset: Subset, p: Price) -> int:
    return min(sum(x[j] for j in subset) for x in demand(i, p))


def over_demand(subset: Subset, p: Price) -> int:
    return sum(mu(i, subset, p) for i in range(3)) - sum(SUPPLY[j] for j in subset)


def inclusion_minimal_maximizer(p: Price) -> Subset:
    values = {x: over_demand(x, p) for x in SUBSETS}
    maximum = max(values.values())
    maximizers = [x for x, value in values.items() if value == maximum]
    minimal = [x for x in maximizers if not any(y < x for y in maximizers)]
    assert len(minimal) == 1, (p, values, maximizers, minimal)
    selected = minimal[0]
    assert all(selected <= x for x in maximizers)
    return selected


def xstar_from_minimal_profile(profile: Sequence[FrozenSet[Bundle]]) -> Subset:
    values = {
        subset: sum(min(sum(x[j] for j in subset) for x in profile[i]) for i in range(3))
        - sum(SUPPLY[j] for j in subset)
        for subset in SUBSETS
    }
    maximum = max(values.values())
    maximizers = [x for x, value in values.items() if value == maximum]
    minimal = [x for x in maximizers if not any(y < x for y in maximizers)]
    assert len(minimal) == 1, (profile, values)
    return minimal[0]


def z_profile_points(p: Price) -> FrozenSet[tuple[int, int]]:
    return frozenset(
        (
            SUPPLY[0] - sum(xs[i][0] for i in range(3)),
            SUPPLY[1] - sum(xs[i][1] for i in range(3)),
        )
        for xs in product(*(demand(i, p) for i in range(3)))
    )


def convex_hull_2d(points: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    pts = sorted(set(points))
    if len(pts) <= 1:
        return tuple(pts)

    def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[int, int]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[tuple[int, int]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    hull = lower[:-1] + upper[:-1]
    # Collinear input gives its two endpoints.
    return tuple(dict.fromkeys(hull))


def check_m_natural_exchange() -> None:
    for i in range(3):
        for x in BUNDLES:
            for y in BUNDLES:
                for j in range(2):
                    if x[j] <= y[j]:
                        continue
                    candidates: list[int | None] = [None] + [k for k in range(2) if x[k] < y[k]]
                    valid = False
                    for k in candidates:
                        xx, yy = list(x), list(y)
                        xx[j] -= 1
                        yy[j] += 1
                        if k is not None:
                            xx[k] += 1
                            yy[k] -= 1
                        x2, y2 = tuple(xx), tuple(yy)
                        if x2 in BUNDLES and y2 in BUNDLES:
                            if valuation(i, x) + valuation(i, y) <= valuation(i, x2) + valuation(i, y2):
                                valid = True
                                break
                    assert valid, (i, x, y, j)


def check_model_and_minimum_witness() -> None:
    check_m_natural_exchange()
    for i in range(3):
        for j in range(2):
            kappa = F(VALUES[i][j], SLOPES[i][j])
            assert kappa >= 0 and F(SLOPES[i][j]) * kappa == VALUES[i][j]

    pbar = (F(5), F(5))
    ds = [demand(i, pbar) for i in range(3)]
    assert (0, 1) in ds[0]
    assert (1, 0) in ds[1]
    assert (0, 0) in ds[2]
    assert tuple(sum(x[j] for x in ((0, 1), (1, 0), (0, 0))) for j in range(2)) == SUPPLY

    # The analytic lower-bound proof uses these exact outside-option thresholds.
    outside_thresholds = (
        (F(10), F(10)),
        (F(10), F(11, 2)),
        (F(0), F(5)),
    )
    assert min(t[1] for t in outside_thresholds) == 5


def check_path_regimes() -> None:
    key = {
        "zero": (F(0), F(0)),
        "quarter": (F(0), F(1, 4)),
        "a": (F(0), F(1, 2)),
        "incoming": (F(1, 2), F(3, 4)),
        "turn": (F(1), F(1)),
        "diag": (F(3), F(3)),
        "minimum": (F(5), F(5)),
    }
    expected = {
        "zero": frozenset({1}),
        "quarter": frozenset({1}),
        "a": frozenset({0, 1}),
        "incoming": frozenset({0, 1}),
        "turn": frozenset({0, 1}),
        "diag": frozenset({0, 1}),
        "minimum": frozenset(),
    }
    for name, p in key.items():
        assert inclusion_minimal_maximizer(p) == expected[name], (name, p)

    assert minimal_bundles(demand(0, key["zero"])) == {(1, 0), (0, 1)}
    assert minimal_bundles(demand(1, key["zero"])) == {(0, 1)}
    assert minimal_bundles(demand(2, key["zero"])) == {(0, 1)}

    # Exact representative grids inside each open segment; the manuscript gives
    # the symbolic utility inequalities for the whole intervals.
    for r in range(1, 20):
        p = (F(0), F(r, 40))
        assert [minimal_bundles(demand(i, p)) for i in range(3)] == [
            frozenset({(1, 0)}), frozenset({(0, 1)}), frozenset({(0, 1)})
        ]
        assert inclusion_minimal_maximizer(p) == frozenset({1})
    for r in range(1, 20):
        t = F(r, 40)
        p = (2 * t, F(1, 2) + t)
        assert [minimal_bundles(demand(i, p)) for i in range(3)] == [
            frozenset({(1, 0)}), frozenset({(1, 0), (0, 1)}), frozenset({(0, 1)})
        ]
        assert inclusion_minimal_maximizer(p) == frozenset({0, 1})
        assert z_profile_points(p) == {(-1, 0), (0, -1)}
    for r in range(1, 40):
        t = F(r, 10)
        if t >= 4:
            break
        p = (1 + t, 1 + t)
        assert [minimal_bundles(demand(i, p)) for i in range(3)] == [
            frozenset({(1, 0), (0, 1)}), frozenset({(1, 0)}), frozenset({(0, 1)})
        ]
        assert inclusion_minimal_maximizer(p) == frozenset({0, 1})


def check_direction_classification() -> None:
    e1, e2, tie = frozenset({(1, 0)}), frozenset({(0, 1)}), frozenset({(1, 0), (0, 1)})
    # p^a: buyer 1=e1, buyer 3=e2, buyer 2 changes according to a1 vs 2a2.
    expected_a = {"<": frozenset({0}), "=": frozenset({0, 1}), ">": frozenset({1})}
    for relation, buyer2 in (("<", e1), ("=", tie), (">", e2)):
        assert xstar_from_minimal_profile((e1, buyer2, e2)) == expected_a[relation]

    # p^circ: a1<a2 gives buyers 1,2 on e1; equality leaves buyer 1 tied;
    # a1>a2 makes buyers 1 and 3 demand e2 regardless of buyer 2's regime.
    assert xstar_from_minimal_profile((e1, e1, e2)) == frozenset({0})
    assert xstar_from_minimal_profile((tie, e1, e2)) == frozenset({0, 1})
    for buyer2 in (e1, tie, e2):
        assert xstar_from_minimal_profile((e2, buyer2, e2)) == frozenset({1})


def check_crossing_certificate() -> None:
    cells = {
        "00": ((F(9, 10), F(1)), (-1, 0)),
        "10": ((F(23, 20), F(11, 10)), (0, -1)),
        "01": ((F(17, 20), F(9, 10)), (0, -1)),
    }
    for name, (p, expected_z) in cells.items():
        assert z_profile_points(p) == {expected_z}, (name, p)

    crossing_points = z_profile_points((F(1), F(1)))
    assert crossing_points == {(-1, 0), (0, -1), (1, -2)}
    assert all(x + y == -1 for x, y in crossing_points)
    assert set(convex_hull_2d(crossing_points)) == {(-1, 0), (1, -2)}

    gradients = [(-1, 0), (0, -1), (0, -2)]
    det = ((gradients[1][0] - gradients[0][0]) * (gradients[2][1] - gradients[0][1])
           - (gradients[1][1] - gradients[0][1]) * (gradients[2][0] - gradients[0][0]))
    assert det == -1


# Optional exact D2 consistency check.  It is not used by the theorem.
def solve_square(rows: Sequence[Sequence[F]], rhs: Sequence[F]) -> tuple[F, ...] | None:
    n = len(rows)
    aug = [list(map(F, rows[i])) + [F(rhs[i])] for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if aug[r][col] != 0), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [v / scale for v in aug[col]]
        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if factor:
                aug[r] = [aug[r][j] - factor * aug[col][j] for j in range(n + 1)]
    return tuple(aug[i][-1] for i in range(n))


def lp_vertices(inequalities: Sequence[tuple[tuple[F, ...], F]], nvars: int) -> FrozenSet[tuple[F, ...]]:
    vertices: set[tuple[F, ...]] = set()
    for active in combinations(range(len(inequalities)), nvars):
        x = solve_square([inequalities[i][0] for i in active], [inequalities[i][1] for i in active])
        if x is not None and all(sum(a[j] * x[j] for j in range(nvars)) >= b for a, b in inequalities):
            vertices.add(x)
    return frozenset(vertices)


def d2_minimizer(delta: F, regime: str) -> tuple[F, FrozenSet[tuple[F, ...]]]:
    e = lambda *xs: tuple(F(x) for x in xs)
    inequalities: list[tuple[tuple[F, ...], F]] = [
        (e(1, 0, 0, 0, 0), F(0)),
        (e(0, 1, 0, 0, 0), F(0)),
    ]
    if regime in {"turn", "outgoing"}:
        inequalities += [(e(0, 0, 1, 0, 0), F(0)), (e(1, 0, 1, 0, 0), F(1)), (e(0, 1, 1, 0, 0), F(1))]
    else:
        inequalities += [(e(0, 0, 1, 0, 0), F(0)), (e(1, 0, 1, 0, 0), F(1))]
    inequalities += [(e(0, 0, 0, 1, 0), F(0)), (e(1, 0, 0, 1, 0), F(1))]
    if regime in {"incoming", "turn"}:
        inequalities.append((e(0, 1, 0, 1, 0), F(1) - delta))
    inequalities += [(e(0, 0, 0, 0, 1), F(0)), (e(0, 1, 0, 0, 1), F(1))]
    vertices = lp_vertices(inequalities, 5)
    objective = lambda x: sum(x)
    value = min(map(objective, vertices))
    return value, frozenset(x for x in vertices if objective(x) == value)


def check_algorithm2_consistency() -> None:
    delta = F(1, 3)  # representative of the symbolic secondary weight gamma*log 2
    value, minimizers = d2_minimizer(delta, "incoming")
    assert value == 2
    assert {(x[0], x[1]) for x in minimizers} == {(F(1), F(1) - delta), (F(1), F(1))}
    for regime in ("turn", "outgoing"):
        value, minimizers = d2_minimizer(delta, regime)
        assert value == 2
        assert {(x[0], x[1]) for x in minimizers} == {(F(1), F(1))}


def main() -> None:
    check_model_and_minimum_witness()
    check_path_regimes()
    check_direction_classification()
    check_crossing_certificate()
    check_algorithm2_consistency()
    print("All exact two-good certificate checks passed.")
    print("M-natural exchange axiom and payment endpoint condition: verified")
    print("path regimes and complete reset-direction classification: verified")
    print("crossing profile points: (-1,0), (0,-1), (1,-2)")
    print("convex-hull extreme points: (-1,0), (1,-2)")
    print("Algorithm-2 dual calculation: optional consistency check passed")
    print("Global equilibrium and convex gluing arguments remain analytic proofs in the manuscript.")


if __name__ == "__main__":
    main()
