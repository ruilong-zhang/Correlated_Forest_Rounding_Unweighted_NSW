# Exact proof functions extracted verbatim; see extraction_manifest.json.
from fractions import Fraction
from flint import arb, fmpq, ctx
def rational(value):
    if isinstance(value,str):
        value=Fraction(value)
    if hasattr(value, 'numerator') and hasattr(value, 'denominator'):
        return fmpq(int(value.numerator), int(value.denominator))
    return fmpq(value)


def exact_endpoint(ball, upper=False):
    endpoint = ball.upper() if upper else ball.lower()
    mantissa, exponent = endpoint.man_exp()
    mantissa, exponent = int(mantissa), int(exponent)
    if exponent >= 0:
        return fmpq(mantissa << exponent)
    return fmpq(mantissa, 1 << -exponent)


def decimal_endpoint(q, digits=24, upper=False):
    q = rational(q)
    scale = 10**digits
    n, d = int(q.p), int(q.q)
    z = -((-n*scale)//d) if upper else (n*scale)//d
    sign = '-' if z < 0 else ''
    z = abs(z)
    return sign+str(z//scale)+'.'+str(z % scale).zfill(digits)


def enclosure(ball, digits=24):
    return [decimal_endpoint(exact_endpoint(ball), digits),
            decimal_endpoint(exact_endpoint(ball, True), digits, True)]


def interval(lo, hi):
    """Ball enclosing the exact rational interval [lo, hi]."""
    lo, hi = rational(lo), rational(hi)
    assert lo <= hi
    return arb((lo+hi)/2, (hi-lo)/2)


def endpoint_geometry(components, p, x0, x1, v0, v1, gap, scale):
    low, diff = arb(0), arb(0)
    constant = True
    for c in components:
        th = c['theta']
        if c['kind'] == 'cap':
            g0 = max(p, min(v0, x0/th), th*gap)
            g1 = max(p, min(v1, x1/th))
        else:
            g0 = max(p, x0, c['A']*v0, th*gap)
            g1 = max(p, x1, c['A']*v1)
        assert 0 < g0 <= g1
        low += arb(c['weight'])*arb(g0/scale).log()
        if g0 != g1:
            constant = False
            diff += arb(c['weight'])*arb(g1/g0).log()
    return low, diff, constant


def maximum_log_ratio(a, v, geo):
    """Certified upper enclosure, including exact endpoint degeneracies.

    Returned ball contains the true maximum whenever the clamp is resolved;
    in the unresolved near-endpoint case its upper endpoint remains valid.
    """
    low, diff, constant = geo
    if a+v == 0:
        return None, 'zero'
    if v == 0 or diff*a >= v:
        return arb(a).log()-low, 'lower'
    if constant or diff*(a+v) <= v:
        return arb(a+v).log()-low-diff, 'upper'
    assert diff > 0 and v > 0
    # The unconstrained stationary maximum is still an upper bound if an
    # endpoint clamp happens to be undecidable at this precision.
    return arb(v).log()-diff.log()-low-1+diff*(a/v), 'interior'
