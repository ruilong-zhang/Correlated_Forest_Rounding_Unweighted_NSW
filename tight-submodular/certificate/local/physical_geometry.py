# Exact arrangements; no discovery cache or floating geometry.
from itertools import combinations
from math import gcd, lcm, comb
from time import monotonic
from fractions import Fraction as Q
from arithmetic import *

P = (1,0,1,0,1)
X0 = (0,1,1,0,0)
V0 = (1,1,1,0,1)
B = (0,0,0,1,1)
X1 = (0,1,1,1,1)
V1 = (1,1,1,1,1)

def enumerate_rays(normals, progress=True):
    n = len(normals)
    wedges = [[None] * n for _ in range(n)]
    coordinates = tuple(combinations(range(5), 2))
    for i, a in enumerate(normals):
        for j in range(i + 1, n):
            b = normals[j]
            wedges[i][j] = tuple(a[u]*b[v] - a[v]*b[u]
                                 for u, v in coordinates)

    rays = set()
    tested = accepted = 0
    start = monotonic()
    # Expand each 4 x 4 determinant into the six products of 2 x 2 minors.
    # Signs alternate with the dropped coordinate, as for a generalized
    # cross product.  Early rejection changes no accepted determinant.
    for i in range(n - 3):
        for j in range(i + 1, n - 2):
            p0,p1,p2,p3,p4,p5,p6,p7,p8,p9 = wedges[i][j]
            for k in range(j + 1, n - 1):
                for q0,q1,q2,q3,q4,q5,q6,q7,q8,q9 in wedges[k][k+1:]:
                    tested += 1
                    v0 = p4*q9-p5*q8+p6*q7+p7*q6-p8*q5+p9*q4
                    sign = (v0 > 0) - (v0 < 0)
                    v1 = -p1*q9+p2*q8-p3*q7-p7*q3+p8*q2-p9*q1
                    if v1:
                        s = 1 if v1 > 0 else -1
                        if sign and s != sign:
                            continue
                        sign = s
                    v2 = p0*q9-p2*q6+p3*q5+p5*q3-p6*q2+p9*q0
                    if v2:
                        s = 1 if v2 > 0 else -1
                        if sign and s != sign:
                            continue
                        sign = s
                    v3 = -p0*q8+p1*q6-p3*q4-p4*q3+p6*q1-p8*q0
                    if v3:
                        s = 1 if v3 > 0 else -1
                        if sign and s != sign:
                            continue
                        sign = s
                    v4 = p0*q7-p1*q5+p2*q4+p4*q2-p5*q1+p7*q0
                    if v4:
                        s = 1 if v4 > 0 else -1
                        if sign and s != sign:
                            continue
                        sign = s
                    if not sign:
                        continue
                    accepted += 1
                    divisor = gcd(v0, v1, v2, v3, v4) * sign
                    rays.add((v0//divisor, v1//divisor, v2//divisor,
                              v3//divisor, v4//divisor))
        if progress and (i % 10 == 0 or i == n - 4):
            print(f'i={i} quadruples={tested} positive={accepted} '
                  f'distinct={len(rays)} seconds={monotonic()-start:.2f}',
                  flush=True)
    assert tested == comb(n, 4), (tested, comb(n, 4))
    return rays, tested, accepted

def primitive(values):
    den = lcm(*(int(v.q) for v in values))
    vec = tuple(int(v*den) for v in values)
    common = gcd(*vec)
    if not common:
        return None
    vec = tuple(v//common for v in vec)
    if next(v for v in vec if v) < 0:
        vec = tuple(-v for v in vec)
    return vec


def joint_normals(components):
    normals = {tuple(int(i==j) for i in range(5)) for j in range(5)}
    for row in components:
        th = rational(row['theta'])
        for upper in (False, True):
            if row['kind'] == 'cap':
                forms = [P, V1 if upper else V0,
                         tuple(fmpq(v)/th for v in (X1 if upper else X0))]
            else:
                A = rational(row['A'])
                forms = [P, X1 if upper else X0,
                         tuple(A*v for v in (V1 if upper else V0))]
            if not upper:
                forms.append(tuple(th*v for v in B))
            for a,b in combinations(forms,2):
                normal = primitive(tuple(rational(x)-rational(y) for x,y in zip(a,b)))
                if normal is not None and (min(normal)<0<max(normal) or sum(v!=0 for v in normal)==1):
                    normals.add(normal)
    return sorted(normals)

def residual_rays(components):
    one,zero=fmpq(1),fmpq(0)
    normals={(1,0,0),(0,1,0),(0,0,1)}
    for c in components:
        th=c['theta']
        for upper in (False,True):
            if c['kind']=='cap':
                forms=[(one,zero,zero),(one,one,one),(zero,1/th,1/th)] if upper else [
                    (one,zero,zero),(one,one,zero),(zero,1/th,zero),(zero,zero,th)]
            else:
                A=c['A']
                forms=[(one,zero,zero),(zero,one,one),(A,A,A)] if upper else [
                    (one,zero,zero),(zero,one,zero),(A,A,zero),(zero,zero,th)]
            for a,b in combinations(forms,2):
                n=primitive(tuple(x-y for x,y in zip(a,b)))
                if n is not None:
                    normals.add(n)
    rays=set()
    for a,b in combinations(sorted(normals),2):
        v=(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
        if not any(v) or min(v)<0<max(v):
            continue
        if max(v)<=0:
            v=tuple(-x for x in v)
        from math import gcd
        common=gcd(*v)
        rays.add(tuple(x//common for x in v))
    return sorted(normals),sorted(rays)

def geometry(base, base_rays, tau):
    if tau is None:
        return base, base_rays, 0
    n = primitive((tau, fmpq(0), tau, fmpq(-1), tau-1))
    normals = sorted(set(base) | {n})
    if tau == 0 or n in base:
        allrays, tested = base_rays, 0
    else:
        allrays, tested, _ = enumerate_rays(normals, progress=False)
    rays = sorted(v for v in allrays if v[3]+v[4] <= tau*(v[0]+v[2]+v[4]))
    assert rays
    return normals, rays, tested

def aq(x):
    return arb(fmpq(str(x)))


def canonical(line):
    line = tuple(Q(x) for x in line)
    lead = next(x for x in line if x)
    return tuple(x/lead for x in line)


def vertices(components):
    # Coordinates are (X/V,p/V).  The domain is the triangle
    # 0<=X,p<=1, X+p>=1.  Redundant lines are harmless.
    lines = {(1,0,0),(1,0,1),(0,1,0),(0,1,1),(1,1,1),(1,-1,0)}
    for item in components:
        if item['kind'] == 'cap':
            theta=Q(item['theta'])
            lines.update([(1,0,theta),(1,-theta,0)])
        else:
            A=Q(item['A'])
            lines.update([(1,0,A),(0,1,A)])
    lines=sorted(set(canonical(row) for row in lines))
    out=set()
    for i,(a,b,c) in enumerate(lines):
        for d,e,f in lines[i+1:]:
            determinant=a*e-b*d
            if not determinant:
                continue
            x,p=(c*e-b*f)/determinant,(a*f-c*d)/determinant
            if 0<=x<=1 and 0<=p<=1 and x+p>=1:
                out.add((x,p))
    return sorted(out),lines


def H(x,p,components):
    result=arb(0)
    for item in components:
        if item['kind']=='cap':
            h=max(p,min(Q(1),x/Q(item['theta'])))
        else:
            h=max(x,p,Q(item['A']))
        assert h>0
        result+=aq(item['weight'])*aq(h).log()
    return result.exp()
