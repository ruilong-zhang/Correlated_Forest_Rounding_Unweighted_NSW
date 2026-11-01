#!/usr/bin/env python3
"""Independent rational/Arb probes for the standalone B implementation.

No bounded cover is run. Checks signed interval arithmetic, log/epsilon
interpolation, the joint score, all credit/price transitions, and boosted
root enclosures (including absent capped roots). Tests supplement the
analytic directed-arithmetic proof; they do not replace global covers.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import argparse
import json
import os
import random
import subprocess
import tempfile
from flint import arb, ctx

ctx.prec=320
SCALE=1<<48


def ar(q):
    q=Q(q);return arb(q.numerator)/q.denominator


def mini(a,b):
    if a<b:return a
    if b<a:return b
    return a.union(b)


def maxi(a,b):
    if a>b:return a
    if b>a:return b
    return a.union(b)


def clip(x,a,b):return mini(maxi(x,a),b)
def dyad(q):
    q=Q(q);return Q(q.numerator*SCALE//q.denominator,SCALE)


z=arb(-1).exp();beta=1-z;ell=arb(2).log()
rs=(1/ell).log()/ell;eps=1+rs-1/ell;c=z-eps;join=beta+rs
a=ar(Q(1,1000));ea=a.exp();cap=(1-(-10*a).exp())/10
lr=ar(Q('1.00032357')).log();pricecap=ar(Q(7,6))*lr
k=1/(4*beta);bcap=ar(Q(43,75))*lr


def logM(D):
    x=ar(Q(9,16))*D
    return (1+mini(x*x/(1+x),ar(Q(5,16))*D)).log()


def logW(T,D):
    if D==0:return maxi(arb(1),T).log()
    r=maxi(arb(0),T/D-beta)
    if r<rs:v=(ell*r).exp()-1
    elif r>rs:v=r-eps
    else:v=((ell*r).exp()-1).union(r-eps)
    return maxi(arb(1),D*(1+v)).log()


def g(s):return mini((1-(-10*s).exp())/10,cap)


def C_ref(t,D,q,mu):
    t,D,q,mu=map(ar,[t,D,q,mu])
    if D==0:return mu*q
    qt=z*t
    if D>=1:onset=beta*D
    elif D<ell:onset=1-c*D
    else:onset=D*(beta-D.log()/ell)
    qs=[q,qt,clip(t-onset,q,qt),clip(t-join*D,q,qt)]
    if D<ea:
        ta=ea-c*D if D<ell*ea else D*(beta+(a-D.log())/ell)
        qs.append(clip(t-ta,q,qt))
    vals=[]
    for qq in qs:
        s=logW(t-qq,D)
        vals.append(mu*qq+mini(logM(D),mini(ar(Q(3,5))*s,s-g(s))))
    value=vals[0]
    for v in vals[1:]:value=mini(value,v)
    return value


def price_ref(w,y,D):
    w,y,D=map(ar,[w,y,D]);s=w.log()
    share=arb(1) if y==0 else mini(arb(1),D/y)
    old=s-mini(logM(w*share),ar(Q(3,5))*s)
    ext=w-1 if w<=1 else g(s)
    return maxi(arb(0),maxi(old-pricecap/2,ar(Q(499,1000))*ext))


def H_ref(t,D,mu):return C_ref(t,D,0,mu)-ar(mu*t)+1+ar(mu).log()


def root_ref(t,D,b):
    floor=k/ar(max(Q(1),D));bar=min(Q(1),1/t) if t else Q(1)
    if H_ref(t,D,bar)<0:return None
    lo=Q(1,100);hi=bar
    if H_ref(t,D,hi)<b:return maxi(floor,ar(hi))
    for _ in range(120):
        mid=(lo+hi)/2;v=H_ref(t,D,mid)-b
        if v<0:lo=mid
        elif v>0:hi=mid
        else:break
    return maxi(floor,ar(lo).union(ar(hi)))


def epsilon_ref(r):
    if r==0:return arb(1)
    r=ar(r);L=(1/r).log1p()
    return r+((1/(r*L)).log()-1)/L


def atomic(path,obj):
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False) as f:
        tmp=Path(f.name);json.dump(obj,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--executable',type=Path,required=True)
    ap.add_argument('--tables',type=Path,required=True)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();assert not args.output.exists()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    work=args.output.parent/(args.output.stem+'_work');assert not work.exists();work.mkdir()
    paths=[args.source,args.executable,args.tables,args.input,Path(__file__)]
    before={str(p.resolve()):sha256(p.read_bytes()).hexdigest() for p in paths}
    tests=[]
    for t in map(Q,['0','.001','.071','.3','.9','1','1.00001','1.2','2']):
        for D in map(Q,['0','.001','.3','.69','.6932','.999','1.0005','1.003','5','20']):
            for mu in map(Q,['.019','.395494177','.625','1']):
                for r in [Q(0),Q(1,10)]:
                    vals=list(map(dyad,[t,D,r*t,mu]));tests.append(('C',vals,C_ref(*vals)))
    for w in map(Q,['.8','1','1.000001','1.0003','1.0008','1.001','1.0015','1.0025','1.5','4','128']):
        for y in map(Q,['0','.1','1','3','64']):
            for D in map(Q,['0','.001','.7','1.2','20']):
                vals=list(map(dyad,[w,y,D]))+[Q(0)];tests.append(('P',vals,price_ref(*vals[:3])))
    fixture=work/'credit_price_inputs.txt'
    fixture.write_text(''.join(kind+' '+' '.join(str(int(v*SCALE)) for v in vals)+'\n' for kind,vals,_ in tests))
    run=subprocess.run([str(args.executable.resolve()),str(args.tables.resolve()),str(args.input.resolve()),'evaluate',str(fixture.resolve())],text=True,capture_output=True,check=True)
    lines=run.stdout.splitlines();assert len(lines)==len(tests)
    maxwidth=Q(0)
    for i,(line,(kind,vals,ref)) in enumerate(zip(lines,tests)):
        token,low,high=line.split();assert token==kind
        lo,hi=Q(int(low),SCALE),Q(int(high),SCALE)
        assert ar(lo)<=ref and ar(hi)>=ref,(i,kind,vals,lo,hi,ref)
        maxwidth=max(maxwidth,hi-lo)
    (work/'credit_price_output.txt').write_text(run.stdout)

    # This harness includes the exact source rather than rewriting any
    # primitive. Its renamed main is unused; it shares the same input data.
    harness=work/'arithmetic_harness.cpp'
    harness.write_text('#define main unused_cover_main\n#include '+json.dumps(str(args.source.resolve()))+'\n#undef main\n'+r'''
int main(int argc,char**argv){
 assert(argc==4);read_tables(argv[1]);read_input(argv[2]);
 std::ifstream f(argv[3]);assert(f);std::string op;i64 a,b,c,d;
 while(f>>op>>a>>b>>c>>d){I v;
  if(op=="J")v=joint_score(PLANAR[0],pt(a),pt(b),pt(c),pt(d));
  else {I x(a,b),y(c,d);
   if(op=="A")v=x+y;else if(op=="S")v=x-y;
   else if(op=="M")v=x*y;else if(op=="D")v=x/y;
   else if(op=="L")v=lg(x);else {assert(op=="E");v=eps_upper(x);}
  }
  std::cout<<op<<" "<<v.l<<" "<<v.h<<"\n";
 }
}
''')
    hexe=work/'arithmetic_harness'
    build=subprocess.run(['g++','-std=c++17','-O3',str(harness),'-o',str(hexe)],capture_output=True,text=True)
    assert build.returncode==0,build.stderr
    rng=random.Random(20261019);primitive=[]
    for _ in range(100):
        nums=[rng.randrange(-32*SCALE,32*SCALE) for _ in range(4)]
        al,ah=sorted(nums[:2]);bl,bh=sorted(nums[2:]);aa,ab,ba,bb=map(lambda x:Q(x,SCALE),[al,ah,bl,bh])
        primitive.extend([('A',[al,ah,bl,bh],(aa+ba,ab+bb)),('S',[al,ah,bl,bh],(aa-bb,ab-ba))])
        products=[x*y for x in (aa,ab) for y in (ba,bb)]
        primitive.append(('M',[al,ah,bl,bh],(min(products),max(products))))
        bl,bh=sorted([rng.randrange(SCALE//16,32*SCALE) for _ in range(2)])
        quotients=[x/Q(y,SCALE) for x in (aa,ab) for y in (bl,bh)]
        primitive.append(('D',[al,ah,bl,bh],(min(quotients),max(quotients))))
    for e in range(-10,11):
        for j in [0,1,511,1024]:
            x=dyad(Q(2)**e*(1+Q(j,1024)));n=int(x*SCALE)
            primitive.append(('L',[n,n,0,0],ar(x).log()))
            # gap_upper bypasses eps_upper for floor/width>=1024.
            # Stay in the reachable epsilon-table domain; a synthetic
            # r=2048 would request the deliberately guarded product 8*r.
            if x<=1024:primitive.append(('E',[n,n,0,0],epsilon_ref(x)))
    cs,ds=map(Q,['.27140539','.72149518']);total=cs+ds
    for _ in range(100):
        t=dyad(Q(rng.randrange(1,2501),1000));D=dyad(Q(rng.randrange(0,21001),1000))
        y=dyad(Q(rng.randrange(1,64001),1000));ss=dyad(min(t,y)*Q(rng.randrange(0,1001),1000))
        mean=ds*y+total*t-ds*ss;width=total*D;floor=max(ds*y,mean-width)
        ref=ar(mean) if width==0 else ar(mean)-ar(width)*epsilon_ref(floor/width)
        primitive.append(('J',[int(v*SCALE) for v in [t,D,y,ss]],ref))
    pinput=work/'arithmetic_inputs.txt'
    pinput.write_text(''.join(op+' '+' '.join(map(str,values))+'\n' for op,values,_ in primitive))
    result=subprocess.run([str(hexe.resolve()),str(args.tables.resolve()),str(args.input.resolve()),str(pinput.resolve())],capture_output=True,text=True,check=True)
    plines=result.stdout.splitlines();assert len(plines)==len(primitive)
    for i,(line,(op,values,ref)) in enumerate(zip(plines,primitive)):
        tag,low,high=line.split();assert tag==op
        lo,hi=Q(int(low),SCALE),Q(int(high),SCALE)
        if isinstance(ref,tuple):assert lo<=ref[0]<=ref[1]<=hi,(op,values,lo,hi,ref)
        else:assert ar(lo)<=ref<=ar(hi),(i,op,values,lo,hi,ref)
    (work/'arithmetic_output.txt').write_text(result.stdout)

    centers=[('0','0'),('.001','0'),('.12','0'),('.13','.000000001'),('.3','.2'),('.5','.5'),
             ('.8','1'),('.999','.001'),('1','.001'),('1.00000001','.001'),('.98','.1'),('1.001','.1'),
             ('1.02','.5'),('1.08','1.08'),('1.0855','1.0855'),('1.10','1.2'),('1.18','1.5'),
             ('1.5','3'),('.5','5'),('.1','20'),('0','20'),('2','5')]
    boxes=[]
    for ts,ds in centers:
        t,D=map(dyad,[ts,ds]);boxes.append((t,t,D,D))
        if t>Q(1,100000) and D>Q(1,100000):
            dt=dyad(Q(1,1000000));boxes.append((t-dt,t+dt,D-dt,D+dt))
    root_replays=excluded=missing=0
    for box in boxes:
        command=[str(args.executable.resolve()),str(args.tables.resolve()),str(args.input.resolve()),'multiplier']+[str(int(x*SCALE)) for x in box]
        output=subprocess.run(command,text=True,capture_output=True,check=True).stdout.split()
        assert output[0]=='MULTIPLIER'
        lo,hi=Q(int(output[1]),SCALE),Q(int(output[2]),SCALE)
        if int(output[3]):excluded+=1;continue
        for t in set(box[:2]):
            for D in set(box[2:]):
                for target in [arb(0),bcap/2,bcap]:
                    ref=root_ref(t,D,target)
                    if ref is None:continue
                    assert ar(lo)<=ref<=ar(hi),(box,t,D,target,lo,hi,ref)
                    root_replays+=1
                    if H_ref(t,D,min(Q(1),1/t) if t else Q(1))<target:missing+=1
    assert root_replays>150 and missing>0
    after={str(p.resolve()):sha256(p.read_bytes()).hexdigest() for p in paths};assert before==after
    report={'status':'PASS_SIMPLE_B_INDEPENDENT_PRIMITIVES','precision_bits':ctx.prec,
            'credit_price_replays':len(tests),'credit_replays':720,'price_replays':275,
            'arithmetic_replays':len(primitive),'root_replays':root_replays,
            'epsilon_probe_domain':'[2^-10,1024], including every block join; larger ratios use gap_upper bypass',
            'root_box_tests':len(boxes),'private_excluded_boxes':excluded,'absent_capped_roots':missing,
            'largest_credit_price_interval':str(maxwidth),'arithmetic_failures':0,
            'pre_sha256':before,'post_sha256':after,
            'harness_sha256':sha256(harness.read_bytes()).hexdigest(),'harness_executable_sha256':sha256(hexe.read_bytes()).hexdigest(),
            'scope':'Independent implementation probes supplement the analytic arithmetic proof. No bounded cover or global approximation is asserted.'}
    atomic(args.output,report)
    print(report['status'],len(tests),'credit/price;',len(primitive),'arithmetic;',root_replays,'roots;',missing,'absent caps; failures 0')


if __name__=='__main__':main()
