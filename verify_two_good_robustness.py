#!/usr/bin/env python3
"""Exact verification for the perturbed two-good obstructions (Lemma H.1).

Companion to verify_two_good_obstruction.py, which fixes the displayed market.
Here the six singleton values and the six payment slopes move over the
admissible class, and every claim of the appendix lemma is rechecked in exact
rational arithmetic:

  conditions (R1)-(R5) of the lemma, evaluated exactly;
  (i)   pbar is the least equilibrium price.  The least element is recomputed
        independently by vertex enumeration over the six assignment-supporting
        polyhedra and compared with the closed form of the lemma;
  (ii)  the auction path: X* on each open segment, the forced incoming ray
        u=(b22,b21), the reset at p^circ, the forced outgoing ray w=(b12,b11),
        and termination at pbar;
  (iii) the three crossing cells, Z(p^circ) a segment, the three limit
        subgradients affinely independent, and the two rays not proportional.

The script also cross-tabulates (R1)-(R5) against the conclusions over a
perturbation battery, which is how the condition list was checked for being
neither too weak nor unnecessarily strong.

Run:  python verify_two_good_robustness.py
"""
from __future__ import annotations

import random
from fractions import Fraction as F
from itertools import combinations, product
from typing import FrozenSet, Sequence

Bundle = tuple[int, int]
Price = tuple[F, F]
Subset = FrozenSet[int]

BUNDLES: tuple[Bundle, ...] = ((0, 0), (1, 0), (0, 1), (1, 1))
SUPPLY: Bundle = (1, 1)
SUBSETS: tuple[Subset, ...] = tuple(
    frozenset(j for j in range(2) if mask & (1 << j)) for mask in range(4)
)

BASE_VALUES = ((F(10), F(10)), (F(10), F(11)), (F(0), F(5)))
BASE_SLOPES = ((F(1), F(1)), (F(1), F(2)), (F(1), F(1)))


class Market:
    """Three unit-demand buyers, two goods in unit supply, linear payments."""

    def __init__(self, values, slopes):
        self.v = tuple(tuple(F(x) for x in row) for row in values)
        self.b = tuple(tuple(F(x) for x in row) for row in slopes)

    def valuation(self, i: int, x: Bundle) -> F:
        if x == (0, 0):
            return F(0)
        return max(self.v[i][j] for j in range(2) if x[j])

    def utility(self, i: int, x: Bundle, p: Price) -> F:
        return self.valuation(i, x) - sum(self.b[i][j] * p[j] * x[j] for j in range(2))

    def demand(self, i: int, p: Price) -> FrozenSet[Bundle]:
        vals = {x: self.utility(i, x, p) for x in BUNDLES}
        best = max(vals.values())
        return frozenset(x for x, u in vals.items() if u == best)

    @staticmethod
    def minimal(ds) -> FrozenSet[Bundle]:
        items = tuple(ds)
        return frozenset(
            x for x in items
            if not any(y != x and all(y[j] <= x[j] for j in range(2)) for y in items)
        )

    def over_demand(self, subset: Subset, p: Price) -> int:
        return (sum(min(sum(x[j] for j in subset) for x in self.demand(i, p))
                    for i in range(3))
                - sum(SUPPLY[j] for j in subset))

    def xstar(self, p: Price) -> Subset:
        vals = {s: self.over_demand(s, p) for s in SUBSETS}
        mx = max(vals.values())
        maximizers = [s for s, val in vals.items() if val == mx]
        minimal = [s for s in maximizers if not any(t < s for t in maximizers)]
        assert len(minimal) == 1, ("no unique minimal maximizer", p, vals)
        return minimal[0]

    def over_demand_max(self, p: Price) -> int:
        return max(self.over_demand(s, p) for s in SUBSETS)

    def z_points(self, p: Price) -> FrozenSet[tuple[int, int]]:
        return frozenset(
            (SUPPLY[0] - sum(xs[i][0] for i in range(3)),
             SUPPLY[1] - sum(xs[i][1] for i in range(3)))
            for xs in product(*(self.demand(i, p) for i in range(3)))
        )

    # ---- closed forms of the lemma --------------------------------------
    def p_a(self) -> Price:
        return (F(0), (self.v[1][1] - self.v[1][0]) / self.b[1][1])

    def p_circ(self) -> Price:
        b = self.b
        det = b[0][0] * b[1][1] - b[0][1] * b[1][0]
        assert det != 0, ("switching lines parallel",)
        r1 = self.v[0][0] - self.v[0][1]
        r2 = self.v[1][0] - self.v[1][1]
        return ((r1 * b[1][1] - r2 * b[0][1]) / det,
                (r1 * b[1][0] - r2 * b[0][0]) / det)

    def p_bar(self) -> Price:
        p2 = self.v[2][1] / self.b[2][1]
        p1 = (self.v[0][0] - self.v[0][1] + self.b[0][1] * p2) / self.b[0][0]
        return (p1, p2)

    # ---- equilibrium price set, computed without the closed form ---------
    def assignment_rows(self, a: int, b_: int):
        """Rows (c1, c2, rhs) meaning c1*p1 + c2*p2 >= rhs for the assignment
        good 1 -> buyer a, good 2 -> buyer b_, third buyer unassigned."""
        c = ({0, 1, 2} - {a, b_}).pop()
        v, B = self.v, self.b
        return [
            (F(1), F(0), F(0)), (F(0), F(1), F(0)),
            (-B[a][0], B[a][1], v[a][1] - v[a][0]), (-B[a][0], F(0), -v[a][0]),
            (B[b_][0], -B[b_][1], v[b_][0] - v[b_][1]), (F(0), -B[b_][1], -v[b_][1]),
            (B[c][0], F(0), v[c][0]), (F(0), B[c][1], v[c][1]),
        ]

    @staticmethod
    def _vertices(rows) -> set[Price]:
        out: set[Price] = set()
        for r1, r2 in combinations(rows, 2):
            det = r1[0] * r2[1] - r1[1] * r2[0]
            if det == 0:
                continue
            p = ((r1[2] * r2[1] - r1[1] * r2[2]) / det,
                 (r1[0] * r2[2] - r1[2] * r2[0]) / det)
            if all(c1 * p[0] + c2 * p[1] >= rhs for c1, c2, rhs in rows):
                out.add(p)
        return out

    def least_equilibrium_price(self) -> Price | None:
        polys = []
        for a, b_ in product(range(3), repeat=2):
            if a == b_:
                continue
            rows = self.assignment_rows(a, b_)
            verts = self._vertices(rows)
            if verts:
                polys.append((rows, verts))
        if not polys:
            return None
        verts = [p for _, vs in polys for p in vs]
        cand = (min(p[0] for p in verts), min(p[1] for p in verts))
        inside = any(all(c1 * cand[0] + c2 * cand[1] >= rhs for c1, c2, rhs in rows)
                     for rows, _ in polys)
        return cand if inside else None

    def xstar_from_profile(self, profile: Sequence[FrozenSet[Bundle]]) -> Subset:
        vals = {s: sum(min(sum(x[j] for j in s) for x in profile[i]) for i in range(3))
                   - sum(SUPPLY[j] for j in s) for s in SUBSETS}
        mx = max(vals.values())
        maximizers = [s for s, val in vals.items() if val == mx]
        minimal = [s for s in maximizers if not any(t < s for t in maximizers)]
        assert len(minimal) == 1, ("no unique minimal maximizer", profile, vals)
        return minimal[0]


E1, E2 = frozenset({(1, 0)}), frozenset({(0, 1)})
TIE = frozenset({(1, 0), (0, 1)})
M = frozenset({0, 1})


def conditions(m: Market) -> str | None:
    """Return None when (R1)-(R5) hold, else the first condition that fails."""
    v, b = m.v, m.b
    if any(v[i][j] < 0 for i in range(3) for j in range(2)):
        return "R1 value negative"
    if any(b[i][j] <= 0 for i in range(3) for j in range(2)):
        return "R1 slope nonpositive"
    if not b[0][1] * b[1][0] < b[1][1] * b[0][0]:
        return "R2 slope ordering"
    try:
        a2 = m.p_a()[1]
        c1, c2 = m.p_circ()
    except (AssertionError, ZeroDivisionError):
        return "R2 switching lines parallel"
    p1, p2 = m.p_bar()
    if not v[1][1] > v[1][0]:
        return "R3 buyer 2 at the zero price"
    if not v[2][0] < v[2][1] - b[2][1] * a2:
        return "R3 buyer 3 on the first segment"
    if not (p1 > 0 and p2 > 0):
        return "R4 pbar not strictly positive"
    if not v[0][0] > b[0][0] * p1:
        return "R4 buyer 1 above the outside option"
    if not v[2][0] < b[2][0] * p1:
        return "R4 buyer 3 rejects good 1"
    if not v[1][0] - b[1][0] * p1 > max(F(0), v[1][1] - b[1][1] * p2):
        return "R4 buyer 2 strictly on good 1"
    if not F(0) < a2 < c2 < p2:
        return "R5 second-coordinate ordering"
    if not F(0) < c1 < p1:
        return "R5 first-coordinate ordering"
    return None


def check_exchange_axiom(m: Market) -> None:
    for i in range(3):
        for x, y in product(BUNDLES, repeat=2):
            for j in range(2):
                if x[j] <= y[j]:
                    continue
                ok = False
                for k in [None] + [k for k in range(2) if x[k] < y[k]]:
                    xx, yy = list(x), list(y)
                    xx[j] -= 1
                    yy[j] += 1
                    if k is not None:
                        xx[k] += 1
                        yy[k] -= 1
                    x2, y2 = tuple(xx), tuple(yy)
                    if x2 in BUNDLES and y2 in BUNDLES and (
                        m.valuation(i, x) + m.valuation(i, y)
                        <= m.valuation(i, x2) + m.valuation(i, y2)
                    ):
                        ok = True
                        break
                assert ok, ("exchange axiom fails", i, x, y, j)


def check_conclusions(m: Market) -> None:
    """Conclusions (i)-(iii) of the lemma."""
    pa, pc, pb = m.p_a(), m.p_circ(), m.p_bar()

    # (i) least equilibrium price, recomputed independently
    least = m.least_equilibrium_price()
    assert least is not None, ("no least equilibrium price",)
    assert least == pb, ("closed form differs from the least element", pb, least)
    assert m.over_demand_max(pb) == 0, ("pbar is not an equilibrium",)

    # (ii) path
    for t in (F(1, 7), F(1, 2), F(6, 7)):
        assert m.xstar((F(0), pa[1] * t)) == frozenset({1}), ("first segment X*",)
    assert m.minimal(m.demand(1, pa)) == TIE, ("buyer 2 not tied at p^a",)
    assert m.minimal(m.demand(0, pa)) == E1, ("buyer 1 not on good 1 at p^a",)
    assert m.xstar(pa) == M, ("X* at p^a",)
    assert m.xstar_from_profile((E1, E1, E2)) == frozenset({0}), ("p^a regime <",)
    assert m.xstar_from_profile((E1, TIE, E2)) == M, ("p^a regime =",)
    assert m.xstar_from_profile((E1, E2, E2)) == frozenset({1}), ("p^a regime >",)

    u = (m.b[1][1], m.b[1][0])
    lam = (pc[0] - pa[0]) / u[0]
    assert lam > 0 and pa[1] + lam * u[1] == pc[1], ("u misses p^circ", u)
    for t in (F(1, 7), F(1, 2), F(6, 7)):
        p = (pa[0] + lam * t * u[0], pa[1] + lam * t * u[1])
        assert [m.minimal(m.demand(i, p)) for i in range(3)] == [E1, TIE, E2], \
            ("incoming regime", p)
        assert m.xstar(p) == M and m.z_points(p) == {(-1, 0), (0, -1)}

    assert m.minimal(m.demand(0, pc)) == TIE, ("buyer 1 not tied at p^circ",)
    assert m.minimal(m.demand(1, pc)) == TIE, ("buyer 2 not tied at p^circ",)
    assert m.minimal(m.demand(2, pc)) == E2, ("buyer 3 not on good 2 at p^circ",)
    assert m.xstar_from_profile((TIE, E1, E2)) == M, ("p^circ regime =",)
    assert m.xstar_from_profile((E1, E1, E2)) == frozenset({0}), ("p^circ regime <",)
    for b2 in (E1, TIE, E2):
        assert m.xstar_from_profile((E2, b2, E2)) == frozenset({1}), ("p^circ regime >",)

    w = (m.b[0][1], m.b[0][0])
    lam2 = (pb[0] - pc[0]) / w[0]
    assert lam2 > 0 and pc[1] + lam2 * w[1] == pb[1], ("w misses pbar", w)
    for t in (F(1, 7), F(1, 2), F(6, 7)):
        p = (pc[0] + lam2 * t * w[0], pc[1] + lam2 * t * w[1])
        assert [m.minimal(m.demand(i, p)) for i in range(3)] == [TIE, E1, E2], \
            ("outgoing regime", p)
        assert m.xstar(p) == M

    # (iii) geometry at the crossing price
    zs = m.z_points(pc)
    assert zs == {(-1, 0), (0, -1), (1, -2)}, ("crossing cells", zs)
    assert all(z1 + z2 == -1 for z1, z2 in zs), ("Z(p^circ) not a segment",)
    g = [(F(-1), F(0)), (F(0), -m.b[0][1] / m.b[0][0]), (F(0), -m.b[1][1] / m.b[1][0])]
    det = ((g[1][0] - g[0][0]) * (g[2][1] - g[0][1])
           - (g[1][1] - g[0][1]) * (g[2][0] - g[0][0]))
    assert det != 0, ("limit subgradients affinely dependent", g)
    assert w[0] * u[1] != w[1] * u[0], ("incoming and outgoing rays proportional",)


def next_event(m: Market, p: Price, a) -> F:
    """Smallest lambda > 0 at which some buyer's demand set changes along a."""
    best = None
    for i in range(3):
        for x in BUNDLES:
            for y in BUNDLES:
                if x >= y:
                    continue
                c0 = m.utility(i, x, p) - m.utility(i, y, p)
                c1 = sum(-m.b[i][j] * a[j] * (x[j] - y[j]) for j in range(2))
                if c1 == 0:
                    continue
                lam = -c0 / c1
                if lam > 0 and (best is None or lam < best):
                    best = lam
    return best if best is not None else F(1)


def check_direction_map(m: Market) -> None:
    """The map from update directions to demand profiles.

    check_conclusions verifies the map from demand profiles to X*.  It never
    verifies which direction produces which profile, so a swapped inequality in
    the proof's direction classification would be invisible to it.  Here a small
    exact step is taken along candidate directions and the resulting profile is
    read off and compared with the classification the proof asserts.
    """
    pa, pc, b = m.p_a(), m.p_circ(), m.b

    def profile_after(p: Price, a):
        eps = next_event(m, p, a) / 2
        q = (p[0] + eps * a[0], p[1] + eps * a[1])
        return tuple(m.minimal(m.demand(i, q)) for i in range(3)), m.xstar(q)

    # at p^a only buyer 2 is tied; the tie survives exactly on u
    u = (b[1][1], b[1][0])
    for a, expect in ((u, TIE), ((u[0] + 1, u[1]), E2), ((u[0], u[1] + 1), E1)):
        sign = b[1][1] * a[1] - b[1][0] * a[0]
        predicted = TIE if sign == 0 else (E1 if sign > 0 else E2)
        assert predicted == expect, ("test direction mislabeled", a)
        prof, xs = profile_after(pa, a)
        assert prof[1] == expect, ("p^a direction map", a, prof)
        assert (xs == M) == (expect is TIE), ("LSC at p^a", a, xs)

    # at p^circ both buyers are tied; thresholds t1 < t2 by (R2)
    t1, t2 = b[0][1] / b[0][0], b[1][1] / b[1][0]
    assert t1 < t2, ("(R2) should order the thresholds", t1, t2)
    for r, expect, xstar in (
        (t1 / 2, (E1, E1, E2), frozenset({0})),
        (t1, (TIE, E1, E2), M),
        ((t1 + t2) / 2, (E2, E1, E2), frozenset({1})),
        (t2, (E2, TIE, E2), frozenset({1})),
        (t2 * 2, (E2, E2, E2), frozenset({1})),
    ):
        prof, xs = profile_after(pc, (r, F(1)))
        assert prof == expect, ("p^circ direction map", r, prof, expect)
        assert xs == xstar, ("p^circ X*", r, xs, xstar)
    w = (b[0][1], b[0][0])
    assert w[0] / w[1] == t1, ("w should sit at the first threshold",)


def perturbed(rng: random.Random, delta: F):
    """Admissible perturbation: schedules stay linear, v31 >= 0."""
    def jitter(x, floor=None):
        y = x + F(rng.randint(-1000, 1000), 1000) * delta
        return max(y, floor) if floor is not None else y
    values = tuple(tuple(jitter(BASE_VALUES[i][j],
                                floor=F(0) if (i, j) == (2, 0) else None)
                         for j in range(2)) for i in range(3))
    slopes = tuple(tuple(jitter(BASE_SLOPES[i][j]) for j in range(2))
                   for i in range(3))
    return values, slopes


def main() -> None:
    base = Market(BASE_VALUES, BASE_SLOPES)
    check_exchange_axiom(base)
    failed = conditions(base)
    assert failed is None, failed
    check_conclusions(base)
    check_direction_map(base)
    print(f"base market: (R1)-(R5) hold, conclusions verified; "
          f"p^a={tuple(str(x) for x in base.p_a())}, "
          f"p^circ={tuple(str(x) for x in base.p_circ())}, "
          f"pbar={tuple(str(x) for x in base.p_bar())}")

    v, b = base.v, base.b
    p1, p2 = base.p_bar()
    print("slacks at the base market: "
          f"R2 {b[1][1] * b[0][0] - b[0][1] * b[1][0]}, "
          f"R3 {v[1][1] - v[1][0]} and {v[2][1] - b[2][1] * base.p_a()[1] - v[2][0]}, "
          f"R4 {v[0][0] - b[0][0] * p1}, {b[2][0] * p1 - v[2][0]}, "
          f"{v[1][0] - b[1][0] * p1 - max(F(0), v[1][1] - b[1][1] * p2)}")

    rng = random.Random(20260902)
    print("\nperturbation battery, (R1)-(R5) against conclusions (i)-(iii):")
    grand = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    for delta in (F(1, 100), F(1, 10), F(1, 2), F(1), F(2)):
        cell = {k: 0 for k in grand}
        for _ in range(600):
            m = Market(*perturbed(rng, delta))
            holds = conditions(m) is None
            try:
                check_conclusions(m)
                if holds:
                    check_direction_map(m)
                concl = True
            except (AssertionError, ZeroDivisionError):
                concl = False
            cell[(holds, concl)] += 1
            grand[(holds, concl)] += 1
        print(f"  delta={str(delta):>6}  hold&true={cell[(True, True)]:4d}  "
              f"hold&FALSE={cell[(True, False)]:4d}  "
              f"fail&true={cell[(False, True)]:4d}  "
              f"fail&false={cell[(False, False)]:4d}")
    assert grand[(True, False)] == 0, ("conditions do not imply the conclusions", grand)
    print(f"\n  totals {grand}")
    print("  (R1)-(R5) imply the conclusions in every draw, as the lemma states.")


if __name__ == "__main__":
    main()
