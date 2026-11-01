#!/usr/bin/env python3
"""Replay the exact local score prerequisites from minimal explicit inputs.

Four Python files total (this entry point plus three extracted proof cores),
parameters.json, and the 4001-node private profile are sufficient. No old
global configuration, score report, ray cache or symmetric archive is read.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
from time import monotonic
import argparse
import json
import os
import sys
import tempfile
from arithmetic import *
from private_profile import Certificate, private_gate
from physical_geometry import residual_rays, joint_normals, enumerate_rays, geometry, vertices, H

HERE=Path(__file__).resolve().parent


def atomic(path,value):
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False) as f:
        tmp=Path(f.name);json.dump(value,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parameters',type=Path,default=HERE/'parameters.json')
    ap.add_argument('--private-profile',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--runtime-generator',type=Path)
    args=ap.parse_args();assert sys.flags.optimize==0
    assert not args.output.exists(),'Do not overwrite a frozen local replay report'
    args.output.parent.mkdir(parents=True,exist_ok=True)
    started=monotonic();ctx.prec=192
    raw=json.loads(args.parameters.read_text());profile=json.loads(args.private_profile.read_text())
    assert raw['components']==profile['components']
    components=[{k:v if k=='kind' else rational(v) for k,v in row.items()} for row in raw['components']]
    assert len(components)==8 and sum(c['weight'] for c in components)==1
    assert all(c['weight']>=0 and 0<c['theta']<=1 for c in components)
    path_checks=[]
    for c in components:
        assert c['kind'] in ('old','cap')
        if c['kind']=='old':
            A,theta=c['A'],c['theta'];assert 0<A<1
            r=int((A/(1-A)).floor());limit=A**(-r)/(r+1)
            assert theta<=limit
            path_checks.append({'A':str(A),'theta':str(theta),'r':r,'limit':str(limit)})
    c0=rational(raw['c0']);assert c0==fmpq(3999999,10000000)
    # The obsolete menu of private affine rows is unnecessary: the final
    # comparison uses c0, the exact mean-profile lower curve, and its
    # analytic log-only tail. This call proves c0 and the harmless A=1 row.
    private_constant=private_gate(components,c0,[fmpq(1)],fmpq(1))
    assert rational(private_constant['strict_log_margin_lower'])>0
    rows=profile['rows'];assert len(rows)==4001
    assert all(Q(row['m'])==Q(i,1000) for i,row in enumerate(rows))
    assert all(0<Q(a['lower_G'])<=Q(b['lower_G']) for a,b in zip(rows,rows[1:]))
    cert=Certificate(raw)
    for i,row in enumerate(rows):
        cert.check(rational(row['m']),rational(row['lower_G']))
        if i%1000==0:print('LOCAL_PRIVATE_PROFILE',i+1,'/',len(rows),'boxes',cert.boxes,flush=True)

    proxies=raw['proxies'];assert len(proxies)==5
    unconditional=[r for r in proxies if Q(r['tau'])==0]
    conditional=[r for r in proxies if Q(r['tau'])>0]
    assert len(unconditional)==1 and len(conditional)==4
    normals,rays=residual_rays(components)
    residual=[];residual_checks=0
    for claim in unconditional:
        assert Q(claim['alpha'])==0
        k,R=rational(claim['k']),rational(claim['R']);upper=None
        for ray in rays:
            p,x,gap=map(rational,ray);scale=rational(max(ray))
            geo=endpoint_geometry(components,p,x,x+gap,p+x,p+x+gap,gap,scale)
            value,_=maximum_log_ratio((p+k*x)/scale,k*gap/scale,geo)
            if value is not None:
                value=exact_endpoint(value,True)
                upper=value if upper is None else max(upper,value)
                residual_checks+=1
        margin=exact_endpoint(arb(R).log())-upper;assert margin>=0
        residual.append({**claim,'log_margin_lower':str(margin)})
    print('LOCAL_UNCONDITIONAL',len(rays),'rays',residual_checks,'checks',flush=True)

    base=joint_normals(components)
    base_rays,base_tested,_=enumerate_rays(base,progress=False)
    base_rays=sorted(base_rays)
    assert all(Q(r['tau'])==4 for r in conditional)
    ns,jrays,tested=geometry(base,base_rays,fmpq(4))
    maxima=[None]*4;conditional_checks=0
    for ray in jrays:
        rr,dd,oo,vv,cc=map(rational,ray)
        p,gap,x0=rr+oo+cc,vv+cc,dd+oo;scale=rational(max(ray))
        geo=endpoint_geometry(components,p,x0,x0+gap,p+dd,p+dd+vv,gap,scale)
        for j,claim in enumerate(conditional):
            alpha,k=map(rational,[claim['alpha'],claim['k']])
            value,_=maximum_log_ratio((p+k*(dd+alpha*oo))/scale,k*(vv+alpha*cc)/scale,geo)
            if value is not None:
                upper=exact_endpoint(value,True)
                maxima[j]=upper if maxima[j] is None else max(maxima[j],upper)
                conditional_checks+=1
    cond=[]
    for claim,upper in zip(conditional,maxima):
        margin=exact_endpoint(arb(rational(claim['R'])).log())-upper
        assert margin>=0
        cond.append({**claim,'log_margin_lower':str(margin)})
    assert len(jrays)==9326 and conditional_checks==37304
    print('LOCAL_CONDITIONAL',len(jrays),'rays',conditional_checks,'checks',flush=True)

    vs,lines=vertices(raw['components']);assert len(vs)==64
    assert len(raw['planar'])==1
    supports=[]
    for row in raw['planar']:
        a,c,d=map(Q,[row['a'],row['c'],row['d']])
        assert a>=0 and c>=0 and d>=0 and 0<a+c+d<=1
        margins=[H(x,p,raw['components'])-arb(rational(a*p+c*x+d)) for x,p in vs]
        assert all(m>0 for m in margins)
        supports.append({**row,'margin_lower':str(min(exact_endpoint(m) for m in margins))})
    handoff=True
    if args.runtime_generator:
        import importlib.util
        spec=importlib.util.spec_from_file_location('lean_runtime_generator',args.runtime_generator)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        actual=[tuple(Q(r[k]) for k in ('tau','alpha','k','R')) for r in proxies]
        assert actual==[tuple(map(Q,r)) for r in module.PROXIES]
        assert [tuple(Q(r[k]) for k in ('a','c','d')) for r in raw['planar']]==[tuple(map(Q,r)) for r in module.PLANAR]
        expected=[(kind,Q(theta),Q(a) if a else None,Q(w)) for kind,theta,a,w in module.COMPONENTS]
        assert [(c['kind'],Q(c['theta']),Q(c['A']) if 'A'in c else None,Q(c['weight'])) for c in raw['components']]==expected
    files=[Path(__file__),HERE/'arithmetic.py',HERE/'private_profile.py',HERE/'physical_geometry.py',args.parameters,args.private_profile]
    if args.runtime_generator:files.append(args.runtime_generator)
    report={'status':'PASS_MINIMAL_LOCAL_CERTIFICATES','global_approximation_claim':False,
            'private_constant':private_constant,'path_checks':path_checks,
            'private_profile':{'points':4001,'interval_boxes':cert.boxes,'maximum_depth':cert.depth},
            'unconditional':{'rays':len(rays),'normals':len(normals),'checks':residual_checks,'rows':residual},
            'conditional':{'base_normals':len(base),'base_rays':len(base_rays),'base_quadruples':base_tested,
                           'normals':len(ns),'rays':len(jrays),'quadruples':tested,'checks':conditional_checks,'rows':cond},
            'joint_planar':{'vertices':len(vs),'lines':len(lines),'supports':supports},
            'runtime_handoff_checked':bool(args.runtime_generator),'seconds':monotonic()-started,
            'input_sha256':{str(p.resolve()):sha256(p.read_bytes()).hexdigest() for p in files},
            'scope':'Exact fixed 8-component private kernel: c0, all 4001 profile nodes, all five affine physical programs, and the single joint support. No historical global claim is used.'}
    atomic(args.output,report)
    print(report['status'],'4001 private nodes;',residual_checks+conditional_checks,'ray checks; 64 planar vertices;',round(report['seconds'],3),'seconds',flush=True)


if __name__=='__main__':main()
