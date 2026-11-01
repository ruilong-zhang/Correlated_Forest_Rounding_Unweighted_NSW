#!/usr/bin/env python3
"""Exact menu compression and B-only NEW redundancy; no scalar cover.

Fractions prove the algebra; Arb certifies the six strict
one-head non-domination examples.  These examples are not approximation
lower bounds.
"""
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse, json, os, sys, tempfile

HERE=Path(__file__).resolve().parent
from flint import arb,ctx
ctx.prec=256
COMPONENT_FILE=HERE/'local/private_profile.json'
CPP=HERE/'verify_B_simple.cpp'

R0,R4,R5,R9,R1=map(Q,['1.00024010','1.00032357','1.02631852','1.19092969','1.23726325'])
k0=Q('.395494177'); c0=Q('.3999999')
cs,ds=map(Q,['.27140539','.72149518']); bs=cs+ds
checks=[]
def ck(name,condition):
    assert condition,name
    checks.append(name)
def ar(q):
    q=Q(q);return arb(q.numerator)/q.denominator
def amax(values):
    values=list(values);out=values[0]
    for v in values[1:]:
        if v>out:out=v
        elif not out>v:out=out.union(v)
    return out

def certify(component_file=COMPONENT_FILE, cpp=CPP):
    da=Q(21,20)*R4*R4
    ck('R4_above_one',R4>1)
    ck('NEW_critical_t_strictly_inside_A',Q(3,40)*R4*R4<Q(3,25))
    ck('NEW_critical_D_inside_proved_A_padding',da<Q('1.05068'))
    ck('NEW_old_eligibility_inclusion',Q(21,20)<4)
    ck('high_t_NEW_head_above_R_squared',Q(40,3)*Q(3,25)>R4*R4)
    ck('high_D_NEW_head_boundary_exact',Q(20,21)*da==R4*R4)
    ck('raw_price_floor_redundant_on_private_domain',c0>0)
    ck('zero_singleton_affine_below_private',k0/R0<c0)

    # A family of lines with increasing slopes and increasing consecutive
    # crossings has exactly this upper-envelope order on the real line.
    pure=[(Q(0),Q(1)),(k0/R0,1/R0),(Q(2,5)/R4,1/R4),(Q(1,2)/R5,1/R5)]
    mixed=[(Q(0),Q(1)),(Q(9,10)/R9,1/R9),(1/R1,1/R1)]
    breakpoints={}
    for name,lines in [('pure',pure),('mixed',mixed)]:
        xs=[]
        for j in range(len(lines)-1):
            sl,cl=lines[j];sr,cr=lines[j+1]
            ck(f'{name}_slope_{j}',sl<sr)
            x=(cl-cr)/(sr-sl);xs.append(x)
            ck(f'{name}_crossing_identity_{j}',cl+sl*x==cr+sr*x)
            ck(f'{name}_positive_crossing_{j}',x>0)
        for j in range(len(xs)-1):ck(f'{name}_ordered_crossing_{j}',xs[j]<xs[j+1])
        for j,x in enumerate(xs):
            v=lines[j][1]+lines[j][0]*x
            ck(f'{name}_no_hidden_line_at_crossing_{j}',all(c+s*x<=v for s,c in lines))
        breakpoints[name]=[str(x) for x in xs]

    components=json.loads(component_file.read_text())['components']
    loggamma=arb(0)
    for c in components:
        loggamma+=ar(c['weight'])*ar(c['theta']).log()
    ck('Gamma_below_one_twentieth',loggamma.exp()<ar(Q(1,20)))
    programs=[('unconditional',Q(0),k0,R0),('d04',Q(0),Q(2,5),R4),
              ('d05',Q(0),Q(1,2),R5),('raw09',Q(1,10),Q(9,10),R9),
              ('raw10',Q(1,10),Q(1),R1)]
    # p,alpha,F,d; none is an eligible NEW row.  HOME is a different head.
    fixtures=[('unconditional',1,5,Q(1,100),Q(1,100)),
              ('d04',1,4,Q(1,10),Q(1,10)),
              ('d05',1,4,Q(2,5),Q(2,5)),
              ('raw09',1,4,Q(11,20),Q(11,20)),
              ('raw10',1,4,Q(7,10),Q(7,10)),
              ('joint',1,Q(1,100),1,1)]
    witnesses=[]
    for winner,p,alpha,F,d in fixtures:
        p,alpha,F,d=map(Q,[p,alpha,F,d])
        ck(winner+'_valid_means',0<=d<=F<=p+d)
        ck(winner+'_NEW_disabled',not(F<=Q(3,40)*p and alpha<=Q(21,20)*p))
        scores={'singleton':ar(p)}
        for name,chi,k,R in programs:
            if name=='unconditional' or alpha<=4*p:
                scores[name]=ar((p+k*(chi*F+(1-chi)*d))/R)
        m=F/alpha
        # An explicit admissible private law supplies an UPPER bound.
        scores['private_upper']=ar(F) if m>=1 else ar(alpha)*(ar(1-m)*loggamma).exp()
        mean=cs*F+ds*(p+d);width=bs*alpha;floor=max(ds*p,mean-width)
        rr=ar(floor/width);ell=(1+1/rr).log()
        epsilon=rr+((1/(rr*ell)).log()-1)/ell
        scores['joint']=ar(mean)-ar(width)*epsilon
        other=amax(v for name,v in scores.items() if name!=winner)
        margin=scores[winner]-other
        ck(winner+'_strictly_above_all_other_original_options',margin>0)
        witnesses.append({'winner':winner,'p':str(p),'alpha':str(alpha),'F':str(F),'d':str(d),
                          'strict_margin':str(margin),'scores':{k:str(v) for k,v in scores.items()}})
    source=cpp.read_text()
    ck('simple_B_has_eight_price_terms','return 2+int(PROX.size()+PLANAR.size())' in source)
    ck('simple_B_omits_NEW_sqrt','new_index' not in source and 'sqroot' not in source)
    ck('simple_B_omits_physical_private_guard','C0*t' not in source)
    ck('simple_B_retains_private_feasibility','if((C0*pt(t.l)).l>SCALE)return true;' in source)
    return {'status':'PASS_EXACT_MENU_COMPRESSION_AND_B_REDUNDANCY',
            'exact_check_count':len(checks),'checks':checks,'D_A':str(da),
            'affine_profile_breakpoints':breakpoints,'strict_nondomination_witnesses':witnesses,
            'source_sha256':sha256(cpp.read_bytes()).hexdigest(),
            'component_sha256':sha256(component_file.read_bytes()).hexdigest(),
            'scope':'Exact same-function affine compression and removal of NEW/private-guard terms from B comparison under the canonical rectangle. Six strict local non-domination examples exclude deletion by pointwise score dominance only; no approximation lower bound and no ablated full cover.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--components',type=Path,default=COMPONENT_FILE)
    parser.add_argument('--source',type=Path,default=CPP)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=certify(args.components,args.source);outpath=args.output
    with tempfile.NamedTemporaryFile('w',dir=outpath.parent,prefix=outpath.name+'.',delete=False) as out:
        tmp=Path(out.name);json.dump(report,out,indent=2);out.write('\n');out.flush();os.fsync(out.fileno())
    os.replace(tmp,outpath)
    print(report['status'],report['exact_check_count'],'exact checks')
    print('D_A =',report['D_A'])
    print('affine breakpoints:',report['affine_profile_breakpoints'])
    for row in report['strict_nondomination_witnesses']:
        print(row['winner'],'strict winner margin',row['strict_margin'])
