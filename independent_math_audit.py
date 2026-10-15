#!/usr/bin/env python3
"""Independent exact audit of the R6 finite-market certificates.

This script was written independently of the two certificate checkers.
It checks only finite/combinatorial and linear-algebra claims.  It does not
replace the convex-analysis proofs (polyhedral gluing, local boundedness, or
closedness of the subdifferential graph).
"""
from fractions import Fraction as Q
from itertools import product

# ---------- Two-good lower-box/path market ----------
B = [(0,0),(1,0),(0,1),(1,1)]
V = [(10,10),(10,11),(0,5)]
S = [(1,1),(1,2),(1,1)]
supply=(1,1)

def val(i,x):
    if x==(0,0): return Q(0)
    return max(Q(V[i][j]) for j in range(2) if x[j])

def util(i,x,p):
    return val(i,x)-sum(Q(S[i][j])*p[j]*x[j] for j in range(2))

def dem(i,p):
    us=[util(i,x,p) for x in B]
    m=max(us)
    return {x for x,u in zip(B,us) if u==m}

def minimal_bundles(D):
    return {x for x in D if not any(y!=x and all(y[j]<=x[j] for j in range(2)) for y in D)}

def mu(i,X,p):
    return min(sum(x[j] for j in X) for x in dem(i,p))

def od(X,p):
    return sum(mu(i,X,p) for i in range(3))-sum(supply[j] for j in X)

def od_table(p):
    sets=[frozenset(),frozenset({0}),frozenset({1}),frozenset({0,1})]
    return {X:od(X,p) for X in sets}

def unique_inclusion_minimal_maximizer(p):
    tab=od_table(p); mx=max(tab.values())
    maxs=[X for X,v in tab.items() if v==mx]
    mins=[X for X in maxs if not any(Y < X for Y in maxs)]
    assert len(mins)==1, (p,tab,maxs,mins)
    return mins[0]

def z_points(p):
    out=set()
    for xs in product(*(dem(i,p) for i in range(3))):
        out.add(tuple(supply[j]-sum(x[j] for x in xs) for j in range(2)))
    return out

# M-natural exchange axiom on the binary domain, with k=None representing 0.
def mnat_check(i):
    for x in B:
        for y in B:
            for j in range(2):
                if x[j] <= y[j]:
                    continue
                candidates=[None]+[k for k in range(2) if x[k]<y[k]]
                ok=False
                for k in candidates:
                    xx=list(x); yy=list(y)
                    xx[j]-=1; yy[j]+=1
                    if k is not None:
                        xx[k]+=1; yy[k]-=1
                    xx=tuple(xx); yy=tuple(yy)
                    if xx in B and yy in B and val(i,x)+val(i,y) <= val(i,xx)+val(i,yy):
                        ok=True; break
                assert ok, (i,x,y,j,candidates)

for i in range(3): mnat_check(i)
assert all(val(i,x)>=0 for i in range(3) for x in B)
# q_ij(kappa)=v_ij: kappa=v/slope, including kappa=0 for buyer 3, good 1.
for i in range(3):
    for j in range(2):
        k=Q(V[i][j],S[i][j])
        assert k>=0 and Q(S[i][j])*k==V[i][j]

# Key demands/minimal demands.
key={
    'zero':(Q(0),Q(0)),
    'quarter':(Q(0),Q(1,4)),
    'a':(Q(0),Q(1,2)),
    'incoming':(Q(1,2),Q(3,4)),
    'turn':(Q(1),Q(1)),
    'diag':(Q(3),Q(3)),
    'min':(Q(5),Q(5)),
}
expected_xstar={
    'zero':frozenset({1}), 'quarter':frozenset({1}),
    'a':frozenset({0,1}), 'incoming':frozenset({0,1}),
    'turn':frozenset({0,1}), 'diag':frozenset({0,1}),
    'min':frozenset(),
}
for name,p in key.items():
    assert unique_inclusion_minimal_maximizer(p)==expected_xstar[name], (name,od_table(p))

assert minimal_bundles(dem(0,key['zero']))=={(1,0),(0,1)}
assert minimal_bundles(dem(1,key['zero']))=={(0,1)}
assert minimal_bundles(dem(2,key['zero']))=={(0,1)}
assert dem(0,key['min']) >= {(1,0),(0,1)}
assert (1,0) in dem(1,key['min'])
assert {(0,0),(0,1)} <= dem(2,key['min'])
assert any(tuple(sum(xs[i][j] for i in range(3)) for j in range(2))==supply
           for xs in product(*(dem(i,key['min']) for i in range(3))))

# Exact path segment profiles, checked at all rational grid points with denominator 40.
for r in range(1,20): # vertical: p2 in (0,1/2)
    p=(Q(0),Q(r,40))
    assert [minimal_bundles(dem(i,p)) for i in range(3)]==[{(1,0)},{(0,1)},{(0,1)}]
    assert unique_inclusion_minimal_maximizer(p)==frozenset({1})
for r in range(1,20): # incoming: (2t,1/2+t), t in (0,1/2)
    t=Q(r,40); p=(2*t,Q(1,2)+t)
    assert minimal_bundles(dem(0,p))=={(1,0)}
    assert minimal_bundles(dem(1,p))=={(1,0),(0,1)}
    assert minimal_bundles(dem(2,p))=={(0,1)}
    assert unique_inclusion_minimal_maximizer(p)==frozenset({0,1})
    assert z_points(p)=={(-1,0),(0,-1)}
for r in range(1,40): # diagonal: (1+t,1+t), t in (0,4)
    t=Q(r,10); 
    if t>=4: break
    p=(1+t,1+t)
    assert minimal_bundles(dem(0,p))=={(1,0),(0,1)}
    assert minimal_bundles(dem(1,p))=={(1,0)}
    assert minimal_bundles(dem(2,p))=={(0,1)}
    assert unique_inclusion_minimal_maximizer(p)==frozenset({0,1})

# Direction classification by all possible strict/tie regimes, not just sample ratios.
# At p^a: b2 chooses 1 if a1<2a2, 2 if >, ties if =.
# Compute X* from the resulting minimal-demand profiles directly.
def xstar_from_profile(profile):
    # profile is list of sets of minimal bundles
    sets=[frozenset(),frozenset({0}),frozenset({1}),frozenset({0,1})]
    vals={X:sum(min(sum(x[j] for j in X) for x in profile[i]) for i in range(3))
             -sum(supply[j] for j in X) for X in sets}
    mx=max(vals.values()); maxs=[X for X,v in vals.items() if v==mx]
    mins=[X for X in maxs if not any(Y<X for Y in maxs)]
    assert len(mins)==1
    return mins[0],vals
base1={(1,0)}; base3={(0,1)}
for rel,b2 in [('<',{(1,0)}),('=',{(1,0),(0,1)}),('>',{(0,1)})]:
    X,tab=xstar_from_profile([base1,b2,base3])
    assert X==({'<':frozenset({0}),'=':frozenset({0,1}),'>':frozenset({1})}[rel]), (rel,tab)
# At turn: if a1<a2 -> b1,b2 good1; if a1=a2 -> b1 ties,b2 good1; if a1>a2 -> b1 good2, b3 good2.
turn_profiles={
    '<':[{(1,0)},{(1,0)},{(0,1)}],
    '=':[{(1,0),(0,1)},{(1,0)},{(0,1)}],
    '>':[{(0,1)},{(1,0)} if False else {(1,0),(0,1)}, {(0,1)}],
}
# For the > case b2 may choose either good depending on a1 vs 2a2; X* is {2} either way.
for prof in ([{(0,1)},{(1,0)},{(0,1)}],[{(0,1)},{(0,1)},{(0,1)}]):
    X,_=xstar_from_profile(prof); assert X==frozenset({1})
X,_=xstar_from_profile(turn_profiles['<']); assert X==frozenset({0})
X,_=xstar_from_profile(turn_profiles['=']); assert X==frozenset({0,1})

# Crossing cells and Z.
cells=[
    ((Q(9,10),Q(1)),(-1,0)),
    ((Q(23,20),Q(11,10)),(0,-1)),
    ((Q(17,20),Q(9,10)),(0,-1)),
]
for p,z in cells:
    assert z_points(p)=={z}
assert z_points(key['turn'])=={(-1,0),(0,-1),(1,-2)}
# The convex hull of these collinear points has endpoints (-1,0),(1,-2), all coordinate sums -1.
assert all(a+b==-1 for a,b in z_points(key['turn']))
# Gradient points are noncollinear for a>0; take a=1.
g=[(-1,0),(0,-1),(0,-2)]
det=(g[1][0]-g[0][0])*(g[2][1]-g[0][1])-(g[1][1]-g[0][1])*(g[2][0]-g[0][0])
assert det==-1

# ---------- Three-good log-affine certificate ----------
K=((1,2,3),(5,2,4),(4,5,1))
ridges={'AC':(15,8,12),'AB':(16,20,30),'BC':(20,25,16)}
expected_ties={
    'AC':({1,2},{0,2},{1}),
    'AB':({1,2},{0},{0,1}),
    'BC':({2},{0,2},{0,1}),
}
for name,y in ridges.items():
    for i in range(3):
        costs=[Q(y[j],K[i][j]) for j in range(3)]
        m=min(costs); arg={j for j,c in enumerate(costs) if c==m}
        assert arg==expected_ties[name][i], (name,i,costs,arg)
# Midpoint single switches.
mid={'A':(Q(31,2),Q(14),Q(21)), 'B':(Q(18),Q(45,2),Q(23)), 'C':(Q(35,2),Q(33,2),Q(14))}
expected_mid={'A':({1,2},{0},{1}), 'B':({2},{0},{0,1}), 'C':({2},{0,2},{1})}
for name,y in mid.items():
    for i in range(3):
        costs=[Q(y[j],K[i][j]) for j in range(3)]
        m=min(costs); arg={j for j,c in enumerate(costs) if c==m}
        assert arg==expected_mid[name][i], (name,i,costs,arg)
# Cell witnesses.
witness={'E':((51,53,58),(2,0,1)), 'A':((31,27,42),(1,0,1)),
         'B':((35,45,46),(2,0,0)), 'C':((36,33,28),(2,2,1))}
for name,(y,profile0) in witness.items():
    got=[]
    for i in range(3):
        costs=[Q(y[j],K[i][j]) for j in range(3)]
        m=min(costs); arg=[j for j,c in enumerate(costs) if c==m]
        assert len(arg)==1; got.append(arg[0])
    assert tuple(got)==profile0,(name,got,profile0)
# Cycle ratios.
from math import prod
lhs=K[0][1]*K[1][2]*K[2][0]
rhs=K[0][2]*K[1][0]*K[2][1]
assert (lhs,rhs)==(32,75)
assert Q(K[1][2],K[0][2])==Q(4,3)
assert Q(K[2][1],K[0][1])*Q(K[1][0],K[2][0])==Q(25,8)

print('INDEPENDENT AUDIT: all finite, exchange-axiom, demand, path, and certificate checks passed.')
print('NOTE: convex-analytic gluing/closed-graph steps remain analytic, not computational.')
