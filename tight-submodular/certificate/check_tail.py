#!/usr/bin/env python3
"""Full-domain boosted singleton tail and every all-count deletion gate."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
from bisect import bisect_right
from time import monotonic
import argparse,json
from scalar import ar,arb,Home,endpoint,decimal


def run(profile):
    start=monotonic();data=json.loads(profile.read_text())
    rows=[(Q(r['m']),Q(r['lower_G'])) for r in data['rows']]
    assert len(rows)==4001 and all(m==Q(i,1000) for i,(m,_) in enumerate(rows))
    assert all(rows[i][1]<=rows[i+1][1] for i in range(4000))
    home=Home(); c0=Q('.3999999'); U0=Q('.6905799'); lam0=Q('.2319048039')
    tmax,Dmax=1/c0,1/rows[0][1]
    assert tmax<3 and Dmax<21 and home.logM(Dmax)<arb(8).log()
    assert 8/(home.z*(1-2*home.z))<119
    floor=endpoint(home.k)/21; k_up=endpoint(home.k,True)
    assert floor>Q(1,100)
    cap_up=endpoint(home.bc,True)
    def tail(m): return max(Q(0),m-1/(8*(m-1))) if m>1 else Q(0)
    def P(m):
        if m>=4: return max(c0*m,rows[-2][1],tail(m))
        i=int(1000*m)
        delayed=rows[max(0,i-1)][1]+(1000*m-i)*(rows[i][1]-rows[max(0,i-1)][1])
        return max(c0*m,delayed,tail(m))
    def private(tl,th,Dl,Dh):
        if Dh==0: return tl
        out=max(c0*tl,Dl*P(tl/Dh))
        if tl>Dh: out=max(out,tl-Dh**2/(8*(tl-Dh)))
        return out
    def upper(tl,th,Dl,Dh):
        geom=home.geometry(tl,Dl)
        allowed=min(Q(1),1/tl) if tl else Q(1)
        lo,hi=floor,allowed
        assert lo<=hi
        for _ in range(24):
            mid=(lo+hi)/2
            test=home.credit(geom,mid)-ar(mid*th)+1+ar(mid).log()
            if test>=ar(cap_up): hi=mid
            else: lo=mid
        mu=min(allowed,max(hi,k_up/max(Q(1),Dl)))
        bound=ar(mu*th)*home.beta+home.logM(Dh)-home.credit(geom,mu)+home.halfK+ar(cap_up)
        return endpoint(bound,True),mu
    stack=[(Q(0),tmax,Q(0),Dmax,0)]
    visits=accept=excluded=depthmax=0;highest=None;witness=None
    transcript=sha256()
    while stack:
        tl,th,Dl,Dh,dep=stack.pop();visits+=1;depthmax=max(depthmax,dep)
        assert visits<100000 and dep<=60,'Unresolved tail cover'
        rec=[str(v) for v in (tl,th,Dl,Dh)]
        pl=private(tl,th,Dl,Dh)
        if pl>1:
            excluded+=1;transcript.update(json.dumps(rec+['private',str(pl)]).encode());continue
        up,mu=upper(tl,th,Dl,Dh)
        if up<U0:
            accept+=1;transcript.update(json.dumps(rec+['tail',str(up),str(mu)]).encode())
            if highest is None or up>highest: highest,witness=up,rec+[str(mu)]
            continue
        tm,Dm=(tl+th)/2,(Dl+Dh)/2
        if (th-tl)*max(Q(1),Dm)>=(Dh-Dl)*max(Q(1),tm):
            stack.extend([(tm,th,Dl,Dh,dep+1),(tl,tm,Dl,Dh,dep+1)])
        else: stack.extend([(tl,th,Dm,Dh,dep+1),(tl,th,Dl,Dm,dep+1)])
    assert highest is not None and visits==2*(accept+excluded)-1
    lam=(arb(4).log()+home.beta.log())/4
    U=(arb(4).log()-3*home.beta.log())/4
    assert U>ar(U0) and lam>ar(lam0)
    deletion={}
    for h,Y in [(2,32),(3,6),(4,2)]:
        gap=(h-1)*(ar(lam0).exp()-1)-(-ar(Q(2,5))*arb(Y).log()+home.halfK).exp()
        assert gap>0
        deletion[str(h)]={'cutoff':Y,'strict_gap':str(gap)}
    zeta=1-(-12*home.L/5).exp()
    high_count={
        'h5':lam-ar(Q(5,4)).log()-ar(Q(6,5))*zeta,
        'derivative':lam-ar(Q(1,5))-ar(Q(6,25))*zeta}
    assert all(v>0 for v in high_count.values())
    return {'status':'PASS_TAIL_AND_ALL_COUNTS','profile_sha256':sha256(profile.read_bytes()).hexdigest(),
      'visited_boxes':visits,'accepted_boxes':accept,'excluded_boxes':excluded,'maximum_depth':depthmax,
      'upper':decimal(highest,30,True),'strict_gap':decimal(U0-highest,30),
      'witness':witness,'transcript_sha256':transcript.hexdigest(),
      'deletions':deletion,'all_count_gates':{k:str(v) for k,v in high_count.items()},
      'elapsed_seconds':monotonic()-start}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('profile',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();r=run(a.profile);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
