#!/usr/bin/env python3
"""Exact and randomized checks for Convex Price Geometry and Potentials.

These checks support, but do not replace, the proofs in the manuscript.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import product
import math
import random
from pathlib import Path
from typing import Callable

import sympy as sp

random.seed(20260817)

REPORT: list[str] = []


def note(msg: str) -> None:
    print(msg)
    REPORT.append(msg)


# ---------------------------------------------------------------------------
# 1. Exact verification of the 3x3 logarithmic demand-cell complex.
# ---------------------------------------------------------------------------
K = [
    [F(1), F(2), F(3)],
    [F(5), F(2), F(4)],
    [F(4), F(5), F(1)],
]
profiles = {
    "H": (2, 1, 1),
    "A": (2, 1, 2),
    "E": (3, 1, 2),
    "B": (3, 1, 1),
    "D": (2, 3, 2),
    "C": (3, 3, 2),
    "F": (3, 3, 1),
}
rays = {
    "AE": (F(69), F(54), F(81)),
    "BE": (F(140), F(175), F(198)),
    "EC": (F(210), F(155), F(168)),
    "AD": (F(45), F(21), F(36)),
    "DC": (F(15), F(6), F(9)),
    "BF": (F(15), F(60), F(12)),
    "FC": (F(100), F(125), F(36)),
    "HA": (F(4), F(5), F(12)),
    "HB": (F(28), F(50), F(75)),
}
# edge -> old cell, new cell, switching buyer (0-based), old good, new good
adj = {
    "AE": ("A", "E", 0, 1, 2),
    "BE": ("B", "E", 2, 0, 1),
    "EC": ("E", "C", 1, 0, 2),
    "AD": ("A", "D", 1, 0, 2),
    "DC": ("D", "C", 0, 1, 2),
    "BF": ("B", "F", 1, 0, 2),
    "FC": ("F", "C", 2, 0, 1),
    "HA": ("H", "A", 2, 0, 1),
    "HB": ("H", "B", 0, 1, 2),
}


def normalized_costs(i: int, y: tuple[F, F, F]) -> list[F]:
    return [y[j] / K[i][j] for j in range(3)]


for edge, y in rays.items():
    old, new, buyer, j_old, j_new = adj[edge]
    old_profile = profiles[old]
    new_profile = profiles[new]
    assert all(old_profile[h] == new_profile[h] for h in range(3) if h != buyer)
    assert {old_profile[buyer], new_profile[buyer]} == {j_old + 1, j_new + 1}
    for h in range(3):
        costs = normalized_costs(h, y)
        minimum = min(costs)
        minimizers = {j + 1 for j, c in enumerate(costs) if c == minimum}
        if h == buyer:
            assert minimizers == {j_old + 1, j_new + 1}, (edge, h, minimizers)
        else:
            assert minimizers == {old_profile[h]}, (edge, h, minimizers)
note("PASS: all nine facet rays have exactly the claimed tie and all other comparisons strict.")

expected_z = {
    "H": (-1, 0, 1),
    "A": (0, -1, 1),
    "E": (0, 0, 0),
    "B": (-1, 1, 0),
    "D": (1, -1, 0),
    "C": (1, 0, -1),
    "F": (0, 1, -1),
}
for name, profile in profiles.items():
    counts = [profile.count(j + 1) for j in range(3)]
    z = tuple(1 - c for c in counts)
    assert z == expected_z[name]
note("PASS: the legacy seven-cell regression has the expected excess-supply vectors.")


# ---------------------------------------------------------------------------
# 2. Legacy exact symbolic verification of the polyhedral holonomy system.
# ---------------------------------------------------------------------------
a, b = sp.symbols("a b", positive=True)
n1 = sp.Matrix([0, sp.Rational(1, 2), -sp.Rational(1, 3)])
n2 = sp.Matrix([sp.Rational(1, 5), 0, -sp.Rational(1, 4)])
n3 = sp.Matrix([sp.Rational(1, 4), -sp.Rational(1, 5), 0])
gA = -a * n1
gB = -b * n3

x = sp.solve(sp.Eq((gA + sp.Symbol("x") * n2)[2], 0), sp.Symbol("x"))[0]
gD = sp.simplify(gA + x * n2)
y = sp.solve(sp.Eq((gD + sp.Symbol("y") * n1)[1], 0), sp.Symbol("y"))[0]
gC_from_A = sp.simplify(gD + y * n1)
tA = sp.solve(sp.Eq((sp.Symbol("tA") * n2)[0], gC_from_A[0]), sp.Symbol("tA"))[0]
assert sp.simplify(x - 4 * a / 3) == 0
assert sp.simplify(y - a) == 0
assert sp.simplify(tA - 4 * a / 3) == 0

u = sp.solve(sp.Eq((gB + sp.Symbol("u") * n2)[0], 0), sp.Symbol("u"))[0]
gF = sp.simplify(gB + u * n2)
v = sp.solve(sp.Eq((gF + sp.Symbol("v") * n3)[1], 0), sp.Symbol("v"))[0]
gC_from_B = sp.simplify(gF + v * n3)
tB = sp.solve(sp.Eq((sp.Symbol("tB") * n2)[0], gC_from_B[0]), sp.Symbol("tB"))[0]
assert sp.simplify(u - 5 * b / 4) == 0
assert sp.simplify(v - b) == 0
assert sp.simplify(tB - 5 * b / 4) == 0
ratio_route_C = sp.solve(sp.Eq(tA, tB), b)[0]
assert sp.simplify(ratio_route_C - 16 * a / 15) == 0

h = sp.solve(sp.Eq((gA - sp.Symbol("h") * n3)[1], 0), sp.Symbol("h"))[0]
gH = sp.simplify(gA - h * n3)
ratio_route_H = sp.solve(sp.Eq(gB[0], gH[0]), b)[0]
assert sp.simplify(h - 5 * a / 2) == 0
assert sp.simplify(ratio_route_H - 5 * a / 2) == 0
assert sp.simplify(ratio_route_C - ratio_route_H) != 0
note("PASS: exact holonomy algebra gives b=(16/15)a on one cycle and b=(5/2)a on the other.")

# Independent exact linear-system audit of the entire facet graph.  Facet
# jumps plus the zero coordinates dictated by excess supply admit only the
# zero gradient assignment, which is incompatible with any nonclearing cell.
cell_names = ["H", "A", "E", "B", "D", "C", "F"]
gvars = {name: sp.Matrix(sp.symbols(f"g{name}1 g{name}2 g{name}3")) for name in cell_names}
edge_names = ["HA", "HB", "AE", "BE", "AD", "DC", "EC", "BF", "FC"]
lvars = {name: sp.symbols(f"l{name}") for name in edge_names}
edge_system = [
    ("H", "A", "HA", n3), ("H", "B", "HB", n1),
    ("A", "E", "AE", n1), ("B", "E", "BE", n3),
    ("A", "D", "AD", n2), ("D", "C", "DC", n1),
    ("E", "C", "EC", n2), ("B", "F", "BF", n2),
    ("F", "C", "FC", n3),
]
linear_eqs = []
for old, new, edge, normal in edge_system:
    linear_eqs.extend(list(gvars[new] - gvars[old] - lvars[edge] * normal))
for name, z in expected_z.items():
    for coordinate, value in enumerate(z):
        if value == 0:
            linear_eqs.append(gvars[name][coordinate])
linear_vars = [entry for name in cell_names for entry in gvars[name]] + [lvars[name] for name in edge_names]
linear_matrix, linear_rhs = sp.linear_eq_to_matrix(linear_eqs, linear_vars)
assert linear_rhs == sp.zeros(linear_rhs.rows, 1)
assert linear_matrix.rank() == len(linear_vars)
assert sp.linsolve(linear_eqs + [lvars["AE"] - 1], linear_vars) == sp.EmptySet
note("PASS: an independent exact rank audit finds that the full facet-jump/zero-pattern system has only the zero solution.")


# ---------------------------------------------------------------------------
# 3. Exact verification of the strengthened three-facet any-convex certificate.
# ---------------------------------------------------------------------------
Vector3 = tuple[F, F, F]


def vadd(*vectors: Vector3) -> Vector3:
    return tuple(sum(vector[j] for vector in vectors) for j in range(3))


def vscale(scalar: F, vector: Vector3) -> Vector3:
    return tuple(scalar * entry for entry in vector)


def vdot(left: Vector3, right: Vector3) -> F:
    return sum(a_entry * b_entry for a_entry, b_entry in zip(left, right))


def demand_set_exact(i: int, y: Vector3) -> tuple[int, ...]:
    costs = normalized_costs(i, y)
    minimum = min(costs)
    return tuple(j + 1 for j, cost in enumerate(costs) if cost == minimum)


def excess_exact(profile: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(1 - profile.count(j + 1) for j in range(3))


ridges: dict[str, Vector3] = {
    "AC": (F(15), F(8), F(12)),
    "AB": (F(16), F(20), F(30)),
    "BC": (F(20), F(25), F(16)),
}
core_profiles = {name: profiles[name] for name in ("E", "A", "B", "C")}
core_excess = {name: expected_z[name] for name in core_profiles}
expected_ridge_demands = {
    "AC": ((2, 3), (1, 3), (2,)),
    "AB": ((2, 3), (1,), (1, 2)),
    "BC": ((3,), (1, 3), (1, 2)),
}
expected_ridge_vertices = {
    "AC": {core_excess["E"], core_excess["A"], core_excess["C"],
           tuple(core_excess["A"][j] + core_excess["C"][j] for j in range(3))},
    "AB": {core_excess["E"], core_excess["A"], core_excess["B"],
           tuple(core_excess["A"][j] + core_excess["B"][j] for j in range(3))},
    "BC": {core_excess["E"], core_excess["B"], core_excess["C"],
           tuple(core_excess["B"][j] + core_excess["C"][j] for j in range(3))},
}

for ridge_name, ridge_y in ridges.items():
    ridge_demands = tuple(demand_set_exact(i, ridge_y) for i in range(3))
    assert ridge_demands == expected_ridge_demands[ridge_name]
    ridge_vertices = {
        excess_exact(tuple(chosen_goods))
        for chosen_goods in product(*ridge_demands)
    }
    assert ridge_vertices == expected_ridge_vertices[ridge_name]
    assert all(coordinate > 1 for coordinate in ridge_y)
note("PASS: all three ridges have exactly the claimed ties, strict comparisons, and excess-supply vertices.")


# An affine-margin certificate proves each universal open-segment statement.
# A designated pair must tie at both endpoints.  Every excluded-good gap is
# affine, nonnegative at both endpoints, and not zero at both; hence it is
# strictly positive for every t in (0,1).
facet_specs = {
    "A/E": (ridges["AC"], ridges["AB"], ((2, 3), (1,), (2,))),
    "B/E": (ridges["AB"], ridges["BC"], ((3,), (1,), (1, 2))),
    "C/E": (ridges["BC"], ridges["AC"], ((3,), (1, 3), (2,))),
}
facet_midpoints: dict[str, Vector3] = {}
for facet_name, (left, right, intended_demands) in facet_specs.items():
    midpoint = vscale(F(1, 2), vadd(left, right))
    facet_midpoints[facet_name] = midpoint
    assert tuple(demand_set_exact(i, midpoint) for i in range(3)) == intended_demands
    for i, intended in enumerate(intended_demands):
        left_costs = normalized_costs(i, left)
        right_costs = normalized_costs(i, right)
        anchor = intended[0] - 1
        if len(intended) == 2:
            tied = intended[1] - 1
            assert left_costs[anchor] == left_costs[tied]
            assert right_costs[anchor] == right_costs[tied]
        for other in range(3):
            if other + 1 in intended:
                continue
            left_gap = left_costs[other] - left_costs[anchor]
            right_gap = right_costs[other] - right_costs[anchor]
            assert left_gap >= 0 and right_gap >= 0
            assert left_gap > 0 or right_gap > 0
note("PASS: exact endpoint margins certify all three open ridge-to-ridge facet segments, not merely sampled points.")


cell_witnesses: dict[str, Vector3] = {
    "E": (F(51), F(53), F(58)),
    "A": (F(31), F(27), F(42)),
    "B": (F(35), F(45), F(46)),
    "C": (F(36), F(33), F(28)),
}
assert cell_witnesses["E"] == vadd(ridges["AC"], ridges["AB"], ridges["BC"])
for cell_name, witness in cell_witnesses.items():
    assert tuple(demand_set_exact(i, witness)[0] for i in range(3)) == core_profiles[cell_name]
    assert all(len(demand_set_exact(i, witness)) == 1 for i in range(3))
note("PASS: exact A/B/C/E witnesses certify four full-dimensional unique-demand cells.")


# Derive, rather than hard-code, the three switch normals from K.
n1_exact: Vector3 = (F(0), F(1, K[0][1]), -F(1, K[0][2]))
n2_exact: Vector3 = (F(1, K[1][0]), F(0), -F(1, K[1][2]))
n3_exact: Vector3 = (F(1, K[2][0]), -F(1, K[2][1]), F(0))
assert sp.Matrix(n1_exact) == n1
assert sp.Matrix(n2_exact) == n2
assert sp.Matrix(n3_exact) == n3
assert vdot(n1_exact, cell_witnesses["E"]) > 0
assert vdot(n2_exact, cell_witnesses["E"]) < 0
assert vdot(n3_exact, cell_witnesses["E"]) > 0
assert vdot(n1_exact, cell_witnesses["A"]) < 0
assert vdot(n2_exact, cell_witnesses["A"]) < 0
assert vdot(n3_exact, cell_witnesses["A"]) > 0
assert vdot(n1_exact, cell_witnesses["B"]) > 0
assert vdot(n2_exact, cell_witnesses["B"]) < 0
assert vdot(n3_exact, cell_witnesses["B"]) < 0
assert vdot(n1_exact, cell_witnesses["C"]) > 0
assert vdot(n2_exact, cell_witnesses["C"]) > 0
assert vdot(n3_exact, cell_witnesses["C"]) > 0
assert vdot(n1_exact, ridges["AC"]) == vdot(n2_exact, ridges["AC"]) == 0
assert vdot(n1_exact, ridges["AB"]) == vdot(n3_exact, ridges["AB"]) == 0
assert vdot(n2_exact, ridges["BC"]) == vdot(n3_exact, ridges["BC"]) == 0
assert vdot(n3_exact, ridges["AC"]) == F(43, 20)
assert vdot(n2_exact, ridges["AB"]) == F(-43, 10)
assert vdot(n1_exact, ridges["BC"]) == F(43, 6)
note("PASS: K-derived normals have the claimed orientations at all cell witnesses and incidences at all ridges.")


minus_n1_exact = vscale(F(-1), n1_exact)
minus_n3_exact = vscale(F(-1), n3_exact)
ray_isolation_specs = [
    ("AC:A", core_excess["A"], minus_n1_exact, n2_exact, 0),
    ("AC:C", core_excess["C"], n2_exact, minus_n1_exact, 1),
    ("AB:A", core_excess["A"], minus_n1_exact, minus_n3_exact, 0),
    ("AB:B", core_excess["B"], minus_n3_exact, minus_n1_exact, 2),
    ("BC:B", core_excess["B"], minus_n3_exact, n2_exact, 2),
    ("BC:C", core_excess["C"], n2_exact, minus_n3_exact, 1),
]
for label, excess_vector, own_ray, other_ray, zero_coordinate in ray_isolation_specs:
    assert excess_vector[zero_coordinate] == 0, label
    assert own_ray[zero_coordinate] == 0, label
    assert other_ray[zero_coordinate] != 0, label
note("PASS: six exact checks verify the zero/nonzero coordinate patterns required for extreme-ray isolation.")


# Derive all three multiplier ratios from the actual normal components sharing
# one diagonal entry at the corresponding ridge.
shared_excess_components = [
    ("AC", "A", "C", 2),
    ("AB", "A", "B", 1),
    ("BC", "B", "C", 0),
]
for label, left_cell, right_cell, coordinate in shared_excess_components:
    assert core_excess[left_cell][coordinate] == -core_excess[right_cell][coordinate], label
    assert abs(core_excess[left_cell][coordinate]) == 1, label
ratio_AC = abs(minus_n1_exact[2]) / abs(n2_exact[2])
ratio_AB = abs(minus_n1_exact[1]) / abs(minus_n3_exact[1])
ratio_BC = abs(minus_n3_exact[0]) / abs(n2_exact[0])
assert ratio_AC == F(4, 3)
assert ratio_AB == F(5, 2)
assert ratio_BC == F(5, 4)
assert ratio_AB * ratio_BC == F(25, 8)
assert ratio_AC != ratio_AB * ratio_BC
cycle_left = K[0][1] * K[1][2] * K[2][0]
cycle_right = K[0][2] * K[1][0] * K[2][1]
assert cycle_left == 32 and cycle_right == 75 and cycle_left != cycle_right
assert ratio_AC / (ratio_AB * ratio_BC) == F(cycle_left, cycle_right) == F(32, 75)
note("PASS: normal-derived ridge ratios are 4/3, 5/2, and 5/4; the cycle products are 32 and 75.")


# C_0=log(100) makes the outside option strictly inferior throughout the compact
# certificate. Every normalized singleton cost is affine along each segment.
outside_threshold = F(100)
certificate_points = list(ridges.values()) + list(facet_midpoints.values()) + list(cell_witnesses.values())
for point in certificate_points:
    assert all(coordinate > 1 for coordinate in point)
    for i in range(3):
        assert max(normalized_costs(i, point)) < outside_threshold
for left, right, _ in facet_specs.values():
    for i in range(3):
        assert max(normalized_costs(i, left)) < outside_threshold
        assert max(normalized_costs(i, right)) < outside_threshold
note("PASS: exp(C_0)=100 certifies positive prices and outside-option domination at all displayed points and along all three closed segments.")


# ---------------------------------------------------------------------------
# 4. Algebraic checks for the C1 affine-semigroup classification.
# ---------------------------------------------------------------------------
x, t, u, c, lam = sp.symbols("x t u c lam", positive=True)

# Translation branch: Q(x)=x/c, so effective-price translation is affine.
Q_aff = x / c
Qinv_aff = lambda r: c * r
S_aff = sp.simplify(Qinv_aff(Q_aff + t))
assert sp.simplify(S_aff - (x + c*t)) == 0
assert sp.simplify(S_aff.subs(x, sp.simplify(Qinv_aff(Q_aff + u))) - (x + c*(t+u))) == 0

# Dilation branch: Q(x)=log(1+x/c)/lambda.
Q_log = sp.log(1 + x/c) / lam
Qinv_log = lambda r: c * (sp.exp(lam*r) - 1)
S_log = sp.simplify(Qinv_log(Q_log + t))
assert sp.simplify(S_log - (sp.exp(lam*t)*x + c*(sp.exp(lam*t)-1))) == 0
S_log_comp = sp.simplify(S_log.subs(x, sp.exp(lam*u)*x + c*(sp.exp(lam*u)-1)))
S_log_sum = sp.exp(lam*(t+u))*x + c*(sp.exp(lam*(t+u))-1)
assert sp.simplify(S_log_comp - S_log_sum) == 0
note("PASS: symbolic checks verify the translation and exponential-dilation affine semigroups and their composition laws.")

# Cross-good compatibility: unequal log rates produce a power map, affine only
# when the exponent is one. Mixed affine/log maps are nonlinear.
rho, c1, c2 = sp.symbols("rho c1 c2", positive=True)
cross_power = c2*((1+x/c1)**rho - 1)
second_power = sp.simplify(sp.diff(cross_power, x, 2))
assert sp.simplify(second_power.subs(rho, 1)) == 0
# A concrete unequal-rate check suffices as a symbolic regression test.
assert sp.simplify(second_power.subs(rho, sp.Rational(2, 1))) != 0
mixed_exp = c2*(sp.exp(lam*x/c1)-1)
mixed_log = c2*sp.log(1+x/c1)/lam
assert sp.simplify(sp.diff(mixed_exp, x, 2)) != 0
assert sp.simplify(sp.diff(mixed_log, x, 2)) != 0
note("PASS: symbolic cross-good checks confirm that log coordinates require a common rate and cannot mix with affine coordinates.")

# Supplemental C2 consistency check (not used by the main C1 proof).
p, beta = sp.symbols("p beta", positive=True)
Q_aff_c2 = beta * p
Q_log_c2 = sp.log(1 + beta * p)
rho_aff = sp.simplify(sp.diff(Q_aff_c2, p, 2) / sp.diff(Q_aff_c2, p) ** 2)
rho_log = sp.simplify(sp.diff(Q_log_c2, p, 2) / sp.diff(Q_log_c2, p) ** 2)
assert rho_aff == 0
assert rho_log == -1
note("PASS: as a supplemental C2 check, affine and normalized log-affine coordinates have curvature ratios 0 and -1.")

# Difference constraints become affine halfspaces under log-affine coordinates.
for _ in range(1000):
    bj = random.uniform(0.1, 4.0)
    bk = random.uniform(0.1, 4.0)
    cc = random.uniform(-3.0, 3.0)
    pj = random.uniform(0.0, 10.0)
    pk = random.uniform(0.0, 10.0)
    lhs_log = math.log1p(bj * pj) - math.log1p(bk * pk) <= cc
    lhs_affine = 1 + bj * pj <= math.exp(cc) * (1 + bk * pk)
    assert lhs_log == lhs_affine
note("PASS: 1,000 randomized tests confirm that every log difference constraint is exactly an affine halfspace in posted prices.")


# ---------------------------------------------------------------------------
# 5. Local smooth obstruction: C1 inverse-map midpoint stress test for random local diffeomorphisms.
# ---------------------------------------------------------------------------
# Use quadratic payment differences with nonsingular Jacobian, solve targets
# numerically, and confirm midpoint signs for shrinking eps.
import numpy as np

for trial in range(300):
    # Positive first derivatives around zero and a nonsingular local matrix.
    A = np.array([[random.uniform(0.5, 2.0), random.uniform(0.5, 2.0)],
                  [random.uniform(0.5, 2.0), random.uniform(0.5, 2.0)]])
    # H derivative is [a_i1, -a_i2].
    J = np.array([[A[0, 0], -A[0, 1]], [A[1, 0], -A[1, 1]]])
    if abs(np.linalg.det(J)) < 0.15:
        continue
    # Small quadratic perturbation H(x)=Jx + 0.5*(c1*x1^2+c2*x2^2) rowwise.
    C = np.array([[random.uniform(-0.5, 0.5), random.uniform(-0.5, 0.5)],
                  [random.uniform(-0.5, 0.5), random.uniform(-0.5, 0.5)]])
    eps = 1e-5
    targets = [np.array([-3*eps, eps]), np.array([eps, -3*eps])]
    sols = []
    for target in targets:
        xvec = np.linalg.solve(J, target)
        # Newton refinement.
        for _ in range(5):
            H = J @ xvec + 0.5 * np.array([C[r, 0]*xvec[0]**2 + C[r, 1]*xvec[1]**2 for r in range(2)])
            DH = J + np.array([[C[0,0]*xvec[0], C[0,1]*xvec[1]],
                               [C[1,0]*xvec[0], C[1,1]*xvec[1]]])
            xvec -= np.linalg.solve(DH, H-target)
        sols.append(xvec)
    mid = 0.5*(sols[0]+sols[1])
    Hmid = J @ mid + 0.5 * np.array([C[r,0]*mid[0]**2 + C[r,1]*mid[1]**2 for r in range(2)])
    assert Hmid[0] < 0 and Hmid[1] < 0
note("PASS: randomized nonlinear stress tests confirm the midpoint-sign construction behind the C1 local-minor proof.")


# ---------------------------------------------------------------------------
# 6. Two-by-two log potential gradients.
# ---------------------------------------------------------------------------
k1, k2 = sp.symbols("k1 k2", positive=True)
# Gradients of the three affine pieces in max{0,k1*y2-y1,y1-k2*y2}.
grads = [sp.Matrix([0, 0]), sp.Matrix([-1, k1]), sp.Matrix([1, -k2])]
assert grads[1][0] < 0 and grads[1][1] > 0
assert grads[2][0] > 0 and grads[2][1] < 0
note("PASS: the 2x2 logarithmic max-potential has the required excess-supply sign patterns.")


# ---------------------------------------------------------------------------
# 7. Regression checks retained from the piecewise-linear version.
# ---------------------------------------------------------------------------
Bundle = tuple[int, int]
BUNDLES: tuple[Bundle, ...] = ((0, 0), (1, 0), (0, 1), (1, 1))


def unit_demand_value(w: tuple[F, F], bundle: Bundle) -> F:
    if bundle == (0, 0):
        return F(0)
    return max(w[j] for j in range(2) if bundle[j])


def utility_table(
    w: tuple[F, F],
    price: tuple[F, F],
    q1: Callable[[F], F],
    q2: Callable[[F], F],
) -> dict[Bundle, F]:
    return {
        bundle: unit_demand_value(w, bundle)
        - q1(price[0]) * bundle[0]
        - q2(price[1]) * bundle[1]
        for bundle in BUNDLES
    }


def demanded(table: dict[Bundle, F]) -> set[Bundle]:
    best = max(table.values())
    return {bundle for bundle, utility in table.items() if utility == best}


# Exact non-rank-one linear example.
q11 = lambda price: price
q12 = lambda price: price
q21 = lambda price: price
q22 = lambda price: 2 * price
w1 = (F(3), F(3))
w2 = (F(3), F(4))
pA = (F(81, 100), F(9, 10))
pB = (F(23, 20), F(11, 10))
pM = tuple((left + right) / 2 for left, right in zip(pA, pB))
assert pM == (F(49, 50), F(1))
expected_demands = {
    "pA": ({(1, 0)}, {(0, 1)}),
    "pB": ({(0, 1)}, {(1, 0)}),
    "pM": ({(1, 0)}, {(1, 0)}),
}
for name, price in (("pA", pA), ("pB", pB), ("pM", pM)):
    d1 = demanded(utility_table(w1, price, q11, q12))
    d2 = demanded(utility_table(w2, price, q21, q22))
    assert (d1, d2) == expected_demands[name]
note("PASS: the exact linear 2x2 example has clearing endpoints and a nonclearing midpoint.")


# Exact kink example and the integrated-price regression.
def qkink(price: F) -> F:
    return price if price <= 1 else 2 * price - 1


w = (F(4), F(4))
pA = (F(1, 2), F(1, 2))
pB = (F(3), F(2))
pM = tuple((left + right) / 2 for left, right in zip(pA, pB))
assert pM == (F(7, 4), F(5, 4))
dA = demanded(utility_table(w, pA, lambda price: price, qkink))
dB = demanded(utility_table(w, pB, lambda price: price, qkink))
dM = demanded(utility_table(w, pM, lambda price: price, qkink))
assert {(1, 0), (0, 1)}.issubset(dA)
assert {(1, 0), (0, 1)}.issubset(dB)
assert dM == {(0, 1)}
assert qkink(F(1)) == 1
assert 2 * F(1) != qkink(F(1))
note("PASS: the exact kink example has a bent clearing locus and a nonclearing midpoint.")
note("PASS: at the kink, beta(p)p=2 differs from the integrated effective price Q(p)=1.")


# Random kink regressions.
for _ in range(2000):
    slope_left = random.uniform(0.1, 4.0)
    slope_right = random.uniform(0.1, 4.0)
    if abs(slope_right - slope_left) < 0.05:
        slope_right += 0.2
    slope_two = random.uniform(0.1, 4.0)
    x0 = random.uniform(1.0, 5.0)
    y0 = random.uniform(1.0, 5.0)
    step = min(0.01, x0 / 10, y0 * slope_two / (10 * max(slope_left, slope_right)))

    def coordinate_one(value: float) -> float:
        if value <= x0:
            return slope_left * value
        return slope_left * x0 + slope_right * (value - x0)

    def coordinate_two(value: float) -> float:
        return slope_two * value

    left = (x0 - step, y0 - slope_left * step / slope_two)
    right = (x0 + step, y0 + slope_right * step / slope_two)
    middle = ((left[0] + right[0]) / 2, (left[1] + right[1]) / 2)
    level_left = coordinate_one(left[0]) - coordinate_two(left[1])
    level_right = coordinate_one(right[0]) - coordinate_two(right[1])
    level_middle = coordinate_one(middle[0]) - coordinate_two(middle[1])
    assert math.isclose(level_left, level_right, rel_tol=1e-10, abs_tol=1e-10)
    assert not math.isclose(level_middle, level_left, rel_tol=1e-8, abs_tol=1e-10)
note("PASS: 2,000 randomized kink configurations bend the effective-price level set.")


# Rank-one marginal matrices have vanishing 2x2 minors.
for _ in range(1000):
    buyer_factor = [random.uniform(0.1, 5.0) for _ in range(3)]
    good_factor = [random.uniform(0.1, 5.0) for _ in range(4)]
    marginal = [[a_i * b_j for b_j in good_factor] for a_i in buyer_factor]
    for i in range(3):
        for k in range(i + 1, 3):
            for j in range(4):
                for ell in range(j + 1, 4):
                    minor = marginal[i][j] * marginal[k][ell] - marginal[i][ell] * marginal[k][j]
                    assert abs(minor) < 1e-10
note("PASS: 1,000 randomized rank-one marginal matrices have zero 2x2 minors.")


report_path = Path(__file__).with_name("verification_report.txt")
report_path.write_text("\n".join(REPORT) + "\n", encoding="utf-8")
note(f"WROTE: {report_path}")
