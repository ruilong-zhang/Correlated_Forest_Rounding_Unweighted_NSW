#!/usr/bin/env python3
"""Exact low/high ordinary and NEW minorants, and the unit completion charge."""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
from time import monotonic
import argparse,json
from scalar import ar,arb,endpoint,decimal


def run(parameters):
    started=monotonic();raw=json.loads(parameters.read_text());cs=raw['components']
    assert len(cs)==8 and sum(Q(c['weight']) for c in cs)==1
    checks={};margins={}
    def ck(name,cond): assert cond,name;checks[name]=True
    def pos(name,value,margin=Q(0)):
        ck(name,value>ar(margin));margins[name]=str(value)
    weights=sum(Q(c['weight'])*Q(c['theta']) for c in cs)
    ck('weighted_protection_floor',weights<Q(37,100))
    caps={Q(c['theta']):Q(c['weight']) for c in cs if c['kind']=='cap'}
    old=[c for c in cs if c['kind']=='old']
    for i,c in enumerate(old):
        a,theta=Q(c['A']),Q(c['theta']);s=int(a/(1-a))
        ck('attenuation_'+str(i),theta*a**s*(s+1)<=1)
    w=2*caps[Q(1)]+caps[Q(1,2)]
    s=3*sum(Q(c['weight'])*(1-Q(c['A']))/Q(c['A']) for c in old)
    ck('unit_charge_w',w<Q(2,3));ck('unit_charge_s',s<Q(1,6))
    slope=Q('.394805930160');off=Q('.000238462746442');A=Q(38,35);cap=Q(19,5)
    def logphysical(y):
        out=arb(0)
        for c in cs:
            th,wt=Q(c['theta']),Q(c['weight'])
            if c['kind']=='cap':value=max(Q(1),min(1+y,y/th))
            else:value=max(Q(1),y,Q(c['A'])*(1+y))
            out+=ar(wt)*ar(value).log()
        return out
    knots={Q(0),cap}
    for c in cs:
        if c['kind']=='cap':
            th=Q(c['theta']);knots.add(th)
            if th<1:knots.add(th/(1-th))
        else:
            a=Q(c['A']);knots.update([1/a-1,a/(1-a)])
    knots=sorted(x for x in knots if 0<=x<=cap)
    ck('low_exact_endpoint_count',len(knots)==12)
    for i,x in enumerate(knots):pos('low_endpoint_'+str(i),logphysical(x)-ar(slope*x-off))
    ck('quartic_cap',cap==Q(7,2)*A)
    ck('quartic_positive_quadratic',Q(1,2)-Q(3)**2/(4*30)==Q(17,40))
    C=47*A**3; d0=Q(15,158); d1=Q(5,23)
    affine=[]
    for k,R in [(Q('.395494177'),Q('1.00024010')),(Q(2,5),Q('1.00032357'))]:
        ck('low_curvature_'+str(k),k*k/(1+k*d0)**2-12*slope*d0*d0/C>Q(7,50))
        def dn(x):return slope*(1+k*x)*(1-4*x**3/C)-k
        lo,hi=Q(0),d0;ck('derivative_bracket_'+str(k),dn(lo)<0<dn(hi))
        for _ in range(96):
            mid=(lo+hi)/2
            if dn(mid)<0:lo=mid
            else:hi=mid
        value=ar(slope*(lo-hi**4/C)-off)-ar(1+k*hi).log()+ar(R).log()
        pos('low_stationary_'+str(k),value,Q(1,400000000))
        ck('wide_curvature_'+str(k),k*k/(1+k*d1)**2-12*slope*d1*d1/C>Q(1,10))
        ck('wide_derivative_'+str(k),dn(d0)>0)
        affine.append({'k':str(k),'bracket':[str(lo),str(hi)],'lower':str(value)})
    # Exact closed-form maximizer of a^j(r+j): j=ceil(a/(1-a)-r), or 0.
    def potential(r):
        if r<=1:return arb(0)
        result=arb(0)
        for c in cs:
            th,wt=Q(c['theta']),Q(c['weight'])
            if c['kind']=='cap':v=max(Q(1),th*r)
            else:
                a=Q(c['A']);q=a/(1-a)-r
                j=max(0,-((-q.numerator)//q.denominator))
                v=max(Q(1),th*a**j*(r+j))
                # Local neighbor tests also catch an erroneous integer ceiling.
                assert v>=th*a**(j+1)*(r+j+1)
                if j:assert v>=th*a**(j-1)*(r+j-1)
            result+=ar(wt)*ar(v).log()
        return result
    def high_kernel(y):
        if 1+y<=A:return ar(1+y).log()+potential(A/(1+y))
        return ar(1+y).log()+potential(A/y)-potential((1+y)/y)
    count=0;minimum=None;left=Q(0);limit=Q(14,3)
    while left<limit:
        right=min(limit,left+Q(1,100))
        value=high_kernel(left)-ar(Q(1,100)+Q(7,20)*right)
        assert value>ar(Q(1,400)),(str(left),str(value))
        low=endpoint(value);minimum=low if minimum is None else min(low,minimum)
        count+=1;left=right
    ck('high_cover_count',count==467)
    pos('high_tail',high_kernel(limit)-ar(Q(1,100)+Q(7,20)*limit),Q(1,400))
    ck('high_marginal_domain',A<Q(30,23)<Q(4,3))
    ck('high_cap',limit==Q(7,2)*Q(4,3))
    C=47*Q(4,3)**3;x=Q(2,23)
    ck('high_decreasing',Q(7,20)*(1+x)<Q(2,5))
    surplus=Q(1,100)+Q(7,20)*(d1-d1**4/C)-(x-x*x/2+x**3/3)
    ck('high_expected_surplus',surplus==Q(678030109,252528518400)>Q(1,500))
    ck('NEW_polynomial',Q(3,10)**2-4*Q(3,5)*Q(1,10)<0)
    pos('NEW_tail',arb(2).log()-ar(Q(27,40)))
    ck('NEW_moment_coefficient',Q(9,10)-Q(3,10)*Q(5,3)==Q(2,5))
    # Wide-region active-program eliminations.
    t,D=Q(1,5),Q(6,5);R4=Q('1.00032357');R5=Q('1.02631852')
    ck('private_inactive',t+Q(37,100)*D<1)
    ck('HOME_inactive',2*t<1)
    ck('raw09_inactive',(Q('1.19092969')-1)/Q('.9')>t)
    ck('raw10_inactive',Q('1.23726325')-1>t)
    ck('joint_inactive',(1-Q('.72149518'))/Q('.99290057')>t)
    ck('point4_enabled',D<4*Q(9,10))
    ck('half_difference_increasing',R4/2-Q(2,5)*R5>0)
    ck('half_dominated',R4-R5+(R4/2-Q(2,5)*R5)*Q(2,9)<0)
    ck('final_head_bounds',1-Q(2,5)*t==Q(23,25) and t/Q(23,25)==d1)
    return {'status':'PASS_COMPLETION_AND_MINORANTS','checks':len(checks),'gates':checks,
      'margins':margins,'unit_charge':{'w':str(w),'s':str(s)},'affine_minima':affine,
      'high_cover':{'intervals':count,'minimum_gap':decimal(minimum),'expected_surplus':str(surplus)},
      'parameters_sha256':sha256(parameters.read_bytes()).hexdigest(),'seconds':monotonic()-started}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parameters',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();r=run(a.parameters);a.output.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['checks'],'gates;',r['high_cover']['intervals'],'high-scale intervals;',round(r['seconds'],3),'seconds')
