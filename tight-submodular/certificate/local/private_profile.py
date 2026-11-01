# The only floating values PROPOSE support slopes; all conclusions are directed.
import math
from arithmetic import *
def rat(x):
    return fmpq(str(x))


class Certificate:
    def __init__(self, data):
        self.components = [(rat(c['theta']), rat(c['weight'])) for c in data['components']]
        assert sum(w for _, w in self.components) == 1
        assert all(0 < t <= 1 and w >= 0 for t, w in self.components)
        self.knots = sorted(set([fmpq(0), fmpq(1)] + [t for t, _ in self.components]))
        self.logs = {t: self.logH(t) for t in self.knots}
        self.boxes = 0
        self.depth = 0
        self.support = self.support_lines()

    def logH(self, x):
        return sum(arb(w) * arb(max(x, t)).log() for t, w in self.components)

    def support_lines(self):
        # Float values only PROPOSE slopes.  Every intercept is rebuilt from
        # exact knot positions and directed lower bounds of every knot value.
        hull = []
        for x in self.knots:
            y = float(self.logs[x])
            while len(hull) >= 2:
                a, ya = hull[-2]
                b, yb = hull[-1]
                if (yb-ya)*(float(x)-b) >= (y-yb)*(b-a):
                    hull.pop()
                else:
                    break
            hull.append((float(x), y))
        lines = []
        for (a, ya), (b, yb) in zip(hull[:-1], hull[1:]):
            slope = fmpq(round((yb-ya)/(b-a)*10**12), 10**12)
            intercept = min(exact_endpoint(self.logs[x]-arb(slope*x)) for x in self.knots)-fmpq(1,10**30)
            # Thus slope*x+intercept <= log H(x) at EVERY kink.
            assert all(arb(slope*x+intercept) <= self.logs[x] for x in self.knots)
            lines.append((slope, intercept))
        return lines

    @staticmethod
    def polynomial_sign(W, m, a, b):
        # J'' has the sign of this cubic on a positive interval.
        power = [W*(m-1), W*(2*m-1), W*(m+1)-(m+2), W-1]
        shifted = []
        for j in range(4):
            shifted.append(sum(power[k]*math.comb(k, j)*a**(k-j)*(b-a)**j
                               for k in range(j, 4)))
        bernstein = [sum(shifted[j]*fmpq(math.comb(k, j), math.comb(3, j))
                         for j in range(k+1)) for k in range(4)]
        if all(v <= 0 for v in bernstein):
            return -1
        if all(v >= 0 for v in bernstein):
            return 1
        return 0

    def fullwidth_branch(self, a, b, m, C, W, target):
        def fun(x):
            low = C if W == 0 else C+arb(W)*x.log()
            return (x+1-arb(m))*low+(arb(m)-x)*(x+1).log()

        def derivative(x):
            if W == 0:
                return C-(1+x).log()+(arb(m)-x)/(1+x)
            return (C+arb(W)*x.log()-(1+x).log()
                    +arb(W)*(1-arb(m)+x)/x+(arb(m)-x)/(1+x))

        # Endpoints are checked separately; concave or monotone intervals
        # therefore need no range enclosure of the function itself.
        fa, fb = fun(arb(a)), fun(arb(b))
        assert fa > target and fb > target, ('endpoint', str(m), str(a), str(b), fa, fb, target)
        stack = [(a, b, 0)]
        while stack:
            u, v, depth = stack.pop()
            self.boxes += 1
            self.depth = max(self.depth, depth)
            curvature = -1 if W == 0 else self.polynomial_sign(W, m, u, v)
            if curvature == -1:
                continue
            box = interval(u, v)
            d = derivative(box)
            if d >= 0 or d <= 0:
                continue
            if fun(box) > target:
                continue
            mid = (u+v)/2
            fm = fun(arb(mid))
            assert fm > target, ('interior target too large', str(m), str(mid), fm, target)
            if curvature == 1:
                # A differentiable convex function is above its tangent.
                # This removes dependency blow-up near a stationary minimum.
                tangent_lower = fm-abs(derivative(arb(mid)))*arb((v-u)/2)
                if tangent_lower > target:
                    continue
            assert depth < 64, ('unresolved branch', str(m), str(u), str(v))
            stack.append((u, mid, depth+1))
            stack.append((mid, v, depth+1))

    def check(self, m, lower):
        assert m >= 0 and lower > 0
        target = arb(lower).log()
        assert self.logH(m) > target
        if m == 0:
            return
        if m <= 1:
            support = max(slope*m+intercept for slope, intercept in self.support)
            assert arb(support) > target, ('kink convex envelope', str(m), str(lower))
        left, right = max(fmpq(0), m-1), m
        cuts = sorted(set([left, right] + [t for t in self.knots if left < t < right]))
        for a, b in zip(cuts[:-1], cuts[1:]):
            mid = (a+b)/2
            W = sum(w for t, w in self.components if t <= mid)
            C = sum(arb(w)*arb(t).log() for t, w in self.components if t > mid)
            self.fullwidth_branch(a, b, m, C, W, target)

def private_gate(components,c0,As,delta):
    assert all(0<A<=1 for A in As)
    knots=sorted({fmpq(0),fmpq(1)}|{c['theta'] for c in components})
    geos=[]
    worst=None
    logthreshold=-arb(c0).log()
    for x in knots:
        low=sum(arb(c['weight'])*arb(max(x,c['theta'])).log() for c in components)
        diff=arb(x+1).log()-low
        constant=(x>=1)
        # At x=1 endpoints1,2, so this is not a zero-diff coin.
        constant=False
        geo=(low,diff,constant)
        value,case=maximum_log_ratio(x,fmpq(1),geo)
        assert value<logthreshold
        u=exact_endpoint(value,True)
        if worst is None or u>worst:worst=u
        geos.append((x,low,diff))
    rows=[]
    for A in As:
        maximum=fmpq(0); witness=None
        for x,low,diff in geos:
            if diff==0:
                q=fmpq(1)
                value=arb(A*(x+q))-(low+diff*q).exp()
            else:
                qball=(arb(A).log()-diff.log()-low)/diff
                if qball<=0:
                    qball=arb(0)
                elif qball>=1:
                    qball=arb(1)
                else:
                    qball=qball.max(arb(0)).min(arb(1))
                value=arb(A)*(arb(x)+qball)-(low+diff*qball).exp()
                q=qball
            u=exact_endpoint(value,True)
            if u>maximum:
                maximum=u
                witness=dict(x=str(x),q=enclosure(arb(q),24))
        # Safe decimal rounding uses exact integers; no float cushion.
        scale=10**10
        rounded=fmpq(-((-int(maximum.p)*scale)//int(maximum.q)),scale)
        assert arb(rounded)-arb(maximum)>=0
        rows.append(dict(A=str(A),D=str(rounded),max_deficit_upper=decimal_endpoint(maximum,24,True),witness=witness))
    a1=next(row for row in rows if row['A']=='1')
    assert rational(a1['D'])<=delta
    logtwo=arb(2).log()
    delta0=1-(1+logtwo.log())/logtwo
    assert delta0<delta
    return dict(c0=str(c0),private_rays=len(knots)+1,
                maximum_log_loss_upper=decimal_endpoint(worst,24,True),
                strict_log_margin_lower=decimal_endpoint(exact_endpoint(logthreshold)-worst,24),
                affine_rows=rows,private_F_minus_delta_B=True,
                home_delta=str(delta),exact_home_deficit=enclosure(delta0,24))
