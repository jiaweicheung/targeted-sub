#!/usr/bin/env python3
"""Exact verification for the three-good log-affine cycle certificate.

The script checks only finite arithmetic and incidence data used by the paper:
  * demand ties and strict inactive comparisons at the three ridges;
  * the three facet midpoints and strict unique-demand cell witnesses;
  * orientation of the clearing cell and active-normal independence;
  * the zero-coordinate/sign compatibility that isolates each scaled
    excess-supply vector on the correct normal ray;
  * the three scale ratios and the cycle mismatch.

It deliberately does NOT assert that the raw excess-supply vectors are
parallel to facet normals.  The paper needs parallelism only after an unknown
positive diagonal scaling.  The convex-analytic endpoint-transport argument
is proved in the manuscript, not delegated to this script.
"""
from __future__ import annotations

from fractions import Fraction as F
from typing import Sequence

K = (
    (F(1), F(2), F(3)),
    (F(5), F(2), F(4)),
    (F(4), F(5), F(1)),
)
RIDGES = {
    "AC": (F(15), F(8), F(12)),
    "AB": (F(16), F(20), F(30)),
    "BC": (F(20), F(25), F(16)),
}
MIDPOINTS = {
    "A": (F(31, 2), F(14), F(21)),
    "B": (F(18), F(45, 2), F(23)),
    "C": (F(35, 2), F(33, 2), F(14)),
}
CELL_WITNESSES = {
    "E": ((F(51), F(53), F(58)), (3, 1, 2)),
    "A": ((F(31), F(27), F(42)), (2, 1, 2)),
    "B": ((F(35), F(45), F(46)), (3, 1, 1)),
    "C": ((F(36), F(33), F(28)), (3, 3, 2)),
}

n1 = (F(0), F(1, 2), F(-1, 3))
n2 = (F(1, 5), F(0), F(-1, 4))
n3 = (F(1, 4), F(-1, 5), F(0))
minus_n1 = tuple(-x for x in n1)
minus_n3 = tuple(-x for x in n3)

zA = (F(0), F(-1), F(1))
zB = (F(-1), F(1), F(0))
zC = (F(1), F(0), F(-1))


def dot(a: Sequence[F], b: Sequence[F]) -> F:
    return sum((x * y for x, y in zip(a, b)), F(0))


def costs(i: int, y: Sequence[F]) -> tuple[F, F, F]:
    return tuple(y[j] / K[i][j] for j in range(3))  # type: ignore[return-value]


def argmin_indices(values: Sequence[F]) -> tuple[int, ...]:
    v = min(values)
    return tuple(j + 1 for j, x in enumerate(values) if x == v)


def midpoint(a: Sequence[F], b: Sequence[F]) -> tuple[F, ...]:
    return tuple((x + y) / 2 for x, y in zip(a, b))


def sign_pattern(x: Sequence[F]) -> tuple[int, ...]:
    return tuple(1 if v > 0 else -1 if v < 0 else 0 for v in x)


def independent(a: Sequence[F], b: Sequence[F]) -> bool:
    # Two nonzero 3-vectors are dependent iff all 2x2 minors vanish.
    return any(a[i] * b[j] - a[j] * b[i] != 0 for i in range(3) for j in range(i + 1, 3))


def check_ridges_and_segments() -> None:
    expected = {
        "AC": ((2, 3), (1, 3), (2,)),
        "AB": ((2, 3), (1,), (1, 2)),
        "BC": ((3,), (1, 3), (1, 2)),
    }
    for name, y in RIDGES.items():
        got = tuple(argmin_indices(costs(i, y)) for i in range(3))
        assert got == expected[name], (name, got, expected[name])
        assert all(c < 100 for i in range(3) for c in costs(i, y))

    # Each open ridge-to-ridge segment lies in one named switching facet.
    segment_data = {
        "A": ("AC", "AB", 0, (2, 3)),
        "B": ("AB", "BC", 2, (1, 2)),
        "C": ("AC", "BC", 1, (1, 3)),
    }
    for facet, (left, right, tying_buyer, tying_goods) in segment_data.items():
        y0, y1 = RIDGES[left], RIDGES[right]
        # The named tie holds at both endpoints, hence on the whole segment.
        for y in (y0, y1):
            c = costs(tying_buyer, y)
            a, b = tying_goods[0] - 1, tying_goods[1] - 1
            assert c[a] == c[b]
        # Every other comparison is either strict at both endpoints or active
        # at only one endpoint; affine gaps are therefore strict in the open segment.
        ym = midpoint(y0, y1)
        expected_mid = {
            "A": ((2, 3), (1,), (2,)),
            "B": ((3,), (1,), (1, 2)),
            "C": ((3,), (1, 3), (2,)),
        }[facet]
        got_mid = tuple(argmin_indices(costs(i, ym)) for i in range(3))
        assert got_mid == expected_mid, (facet, got_mid)
        assert MIDPOINTS[facet] == ym


def check_cell_witnesses() -> None:
    for name, (y, expected_choices) in CELL_WITNESSES.items():
        got = tuple(argmin_indices(costs(i, y)) for i in range(3))
        assert all(len(s) == 1 for s in got), (name, got)
        got_choices = tuple(s[0] for s in got)
        assert got_choices == expected_choices, (name, got_choices, expected_choices)
        assert all(coord > 1 for coord in y)
        assert all(c < 100 for i in range(3) for c in costs(i, y))


def check_orientations_and_ray_isolation() -> None:
    yE = CELL_WITNESSES["E"][0]
    assert dot(n1, yE) > 0
    assert dot(n2, yE) < 0
    assert dot(n3, yE) > 0

    # Positive diagonal scaling preserves zero coordinates and signs.  These
    # patterns are exactly compatible with the target outward normal rays.
    assert sign_pattern(zA) == sign_pattern(minus_n1)
    assert sign_pattern(zB) == sign_pattern(minus_n3)
    assert sign_pattern(zC) == sign_pattern(n2)

    active = {
        "AC": (minus_n1, n2),
        "AB": (minus_n1, minus_n3),
        "BC": (minus_n3, n2),
    }
    for name, (r1, r2) in active.items():
        assert independent(r1, r2), name

    # The zero coordinate of each scaled z-vector eliminates one generator
    # of the two-ray normal cone, exactly as used in Appendix E.
    isolation = [
        ("AC-A", zA, 0, minus_n1, n2),
        ("AC-C", zC, 1, n2, minus_n1),
        ("AB-A", zA, 0, minus_n1, minus_n3),
        ("AB-B", zB, 2, minus_n3, minus_n1),
        ("BC-B", zB, 2, minus_n3, n2),
        ("BC-C", zC, 1, n2, minus_n3),
    ]
    for name, z, coordinate, desired, excluded in isolation:
        assert z[coordinate] == 0
        assert desired[coordinate] == 0
        assert excluded[coordinate] != 0
        assert sign_pattern(z) == sign_pattern(desired), name


def check_cycle_arithmetic() -> None:
    direct = K[1][2] / K[0][2]       # K_23 / K_13 = 4/3
    via_ab = K[2][1] / K[0][1]       # K_32 / K_12 = 5/2
    via_bc = K[1][0] / K[2][0]       # K_21 / K_31 = 5/4
    assert direct == F(4, 3)
    assert via_ab == F(5, 2)
    assert via_bc == F(5, 4)
    assert direct != via_ab * via_bc

    chi = (K[0][1] * K[1][2] * K[2][0]) / (K[0][2] * K[1][0] * K[2][1])
    assert chi == F(32, 75)
    assert chi != 1


def main() -> None:
    check_ridges_and_segments()
    check_cell_witnesses()
    check_orientations_and_ray_isolation()
    check_cycle_arithmetic()
    print("All exact log-affine certificate checks passed.")
    print("ridge ties, open facet segments, and strict cell witnesses: verified")
    print("normal orientations, sign/zero ray isolation: verified")
    print("direct ratio 4/3; indirect ratio 25/8; chi(K)=32/75")
    print("Convex endpoint transport remains an analytic proof in the manuscript.")


if __name__ == "__main__":
    main()
