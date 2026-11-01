#!/usr/bin/env python3
"""Compact global HOME endpoint: full directed base cover and analytic gates."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse, json, sys, time

from flint import arb, ctx
ctx.prec=192

def ar(x):
    x=Q(x);return arb(x.numerator)/x.denominator

def run(profile,output):
    start=time.monotonic()
    data=json.loads(profile.read_text())
    rows=[(Q(x['m']),Q(x['lower_G'])) for x in data['rows']]
    assert len(rows)==4001 and all(rows[i][0]==Q(i,1000) for i in range(4001))
    assert all(rows[i][1]<=rows[i+1][1] for i in range(4000))
    z=arb(-1).exp();beta=1-z;ell=arb(2).log()
    rs=(1/ell).log()/ell;eps=1+rs-1/ell;c=z-eps
    h=beta+rs
    def H(s):
        r=(s-beta).max(arb(0))
        if r < rs:
            return (ell*r).exp()
        if r > rs:
            return 1+r-eps
        # At a branch straddle take the hull, not a guessed branch.
        return (ell*r).exp().union(1+r-eps)
    visits=accepted=maxdepth=0
    smallest=None;witness=None
    for i in range(4000):
        a,pa=rows[i];b=rows[i+1][0]
        assert (b-pa)*(pa+Q(9,16)) < Q(81,256),(i,'M quadratic')
        assert b-pa < Q(5,16),(i,'M linear')
        if b<=pa:
            continue
        upper_log_t=ar(b/pa).log()
        stack=[(Q(0),Q(3,8),0)]
        while stack:
            ql,qh,depth=stack.pop();visits+=1;maxdepth=max(maxdepth,depth)
            r=(H(ar(a*(1-qh)))/ar(pa)).max(arb(1)).log()
            smooth=r+( (-10*r).exp()-1)/10
            value=ar(ql)-upper_log_t+(ar(Q(3,5))*r).min(smooth)
            if value>ar(Q(1,1000)):
                accepted+=1
                if smallest is None or value.lower()<smallest:
                    smallest=value.lower();witness=[i,str(ql),str(qh)]
            else:
                assert depth<20,(i,str(ql),str(qh),str(value),'unresolved')
                mid=(ql+qh)/2
                stack.extend([(ql,mid,depth+1),(mid,qh,depth+1)])
        if i and i%1000==0:
            print('HOME',i,'visited',visits,'seconds',round(time.monotonic()-start,2),flush=True)
    # New analytic m>=4 tail for f10.  Below are exact rational gates.
    assert Q(96,95)*(1-Q(133,144)**10)<1
    x=Q(11,133)
    assert x-x*x/2+x*x*x/3<Q(2,25)
    # (7/2)(7/25)^2(8m(m-1)-1)-(m+2/7)^2, shifted at m=4.
    A=Q(7,2)*Q(7,25)**2
    coeff=[A*95-Q(30,7)**2, A*56-Q(60,7), A*8-1]
    assert min(coeff)>0
    # Unit charge and active-M barrier need only rational gates.
    w=2*Q('.166466453839')+Q('.274867781350')
    s=3*sum(Q(wt)*(1-Q(at))/Q(at) for at,wt in [
        ('.75','.048883852144'),('.80','.044226298717'),
        ('.90','.032834490836'),('.925','.021914693827'),('.999','.254146426697')])
    assert w<Q(2,3) and s<Q(1,6)
    cap=Q(43,75)*Q('.00032357')
    c0=Q('0.3999999')
    assert cap<Q(1,1000) and Q(1000,999)/c0<3
    assert Q(2,5)-cap>Q(1,3) and 1-3*cap>Q(1,3)
    # Smooth branch: its value is >=c_B when r>=.06 and m>=4.
    r0=ar(Q(3,50))
    assert r0+( (-10*r0).exp()-1)/10 > ar(Q(1,95))+ar(cap)
    assert (-10*r0).exp()-ar(cap)>ar(Q(1,3))
    r1=ar(Q(1,100))
    assert r1+((-10*r1).exp()-1)/10>ar(cap)
    assert (-10*r1).exp()-ar(cap)>ar(Q(1,3))
    result={'status':'PASS_COMPACT_HOME_BASE_ENDPOINT',
      'profile_sha256':sha256(profile.read_bytes()).hexdigest(),
      'private_mean_domain':['0','4'],'q_domain':['0','3/8'],
      'profile_intervals':4000,'visited_q_boxes':visits,'accepted_q_boxes':accepted,
      'home_candidate_reserve':'1/1000',
      'maximum_depth':maxdepth,'smallest_box_margin':str(smallest),'witness':witness,
      'tail_shifted_polynomial_coefficients':[str(v) for v in coeff],
      'elapsed_seconds':time.monotonic()-start}
    result['analytic_gates']=constant_gates(profile)
    result['status']='PASS_COMPACT_GLOBAL_HOME'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def constant_gates(profile_path):
    checks,margins={},{}
    def check(name,condition):
        assert condition,name
        checks[name]=True
    def positive(name,value):
        assert value>0,(name,str(value))
        checks[name]=True;margins[name]=str(value)
    z = arb(-1).exp()
    beta = 1-z
    k = 1/(4*beta)
    ell = arb(2).log()
    rs = (1/ell).log()/ell
    eps = 1+rs-1/ell
    c = z-eps
    h = beta+rs
    c0 = ar(".3999999")
    R4 = ar("1.00032357")
    bc = 43*R4.log()/75
    a0 = ar(".001")
    b0 = (1-(-10*a0).exp())/10
    K = 7*R4.log()/6
    check("compact_q_domain_contains_true_domain", 0 < z < ar(Q(3, 8)))
    check("HOME_c_bracket", ar(Q(7, 25)) < c < ar(Q(2, 7)))
    check("tail_is_below_exponential_floor_scale", ar(Q(24, 95)) < ell)
    check("tail_HOME_argument_beyond_affine_join", 4*beta > h)
    check("floor_admissibility_uses_only_c0", k < c0)
    check("global_floor_above_five_candidate_cutoff", k/21 > ar(Q(1, 100)))
    profile = json.loads(profile_path.read_text())
    check("private_floor_bounds_D", Q(profile["rows"][0]["lower_G"]) > Q(1, 21))
    check("tail_smooth_lower_coefficient",
          5-Q(50, 3)*Q(2, 25) >= Q(7, 2))
    check("tail_smooth_implies_linear_HOME",
          5*Q(2, 25) == Q(2, 5) < Q(3, 5))
    check("tail_M_linear_base", Q(1, 24) < Q(5, 16))
    check("tail_M_quadratic_base", Q(5, 24) < Q(1, 4) < Q(81, 256))

    # Active-M elasticity: for x<=5/4, x^3+2x^2<=65x/16.
    check("active_M_crossover", Q(5, 4)**2/(1+Q(5, 4)) == Q(5, 9)*Q(5, 4))
    check("active_M_elasticity_majorant",
          Q(5, 4)**2 + 2*Q(5, 4) == Q(65, 16))
    check("active_M_elasticity_below_two_thirds", Q(5, 4)/16-2 < 0)
    check("linear_M_subcap_scale", 1/(c0*(1-bc)) < 3)
    check("M_subcap_slope_below_t_one", 1-3*bc > ar(Q(1, 3)))
    check("linear_HOME_subcap_slope", ar(Q(2, 5))-bc > ar(Q(1, 3)))
    for label, r, ceiling in [
        ("small_r", ar(Q(1, 100)), bc),
        ("tail_r", ar(Q(3, 50)), bc+ar(Q(1, 95))),
    ]:
        f = r + ((-10*r).exp()-1)/10
        positive("smooth_subcap_height_"+label, f-ceiling)
        positive("smooth_subcap_derivative_"+label, (-10*r).exp()-bc-ar(Q(1, 3)))
    check("compact_start_above_boost_cap", bc < ar(Q(1, 1000)))
    check("zero_D_near_one", 1-bc > ar(Q(1, 3)))

    # New exact small-scale tail credit.
    A, dmax, rate = ar(Q(81, 256)), ar(Q(11, 1000)), ar(Q(27, 100))
    positive("tail_W_linear_gain", ar(Q(7, 25))-A*dmax/k-rate)
    denom = 1+rate*dmax
    positive("tail_smooth_credit_coefficient",
             (5-ar(Q(50, 3))*rate*dmax)*(rate/denom)**2-A)
    positive("tail_linear_credit_coefficient",
             ar(Q(3, 5))*rate/denom-A*dmax)
    positive("localized_five_candidate_derivative",
             ar(Q(1, 100))-(1-(-10*a0).exp()))
    positive("localized_inactive_shortcut",
             ar(Q(2,5))*ar(Q(401,400)).log()-b0)
    positive("half_price_net_nonnegative",
             ar(Q(501, 1000))*(1-5*K)-ar(Q(1, 2)))
    check("half_price_cap_below_join", K < a0)
    check("localized_max_credit_decrease", b0-ar(Q(2, 5))*a0 < ar(Q(3, 5000)))

    return {'count':len(checks),'checks':checks,'margins':margins}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('profile',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();run(args.profile,args.output)
