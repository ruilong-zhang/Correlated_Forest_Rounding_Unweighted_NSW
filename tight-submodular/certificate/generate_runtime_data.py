#!/usr/bin/env python3
"""Generate only the directed tables and exact rational input used by B.

Runtime dependencies: Python standard library, python-flint, and the
certified 4001-node private profile JSON supplied by --private-profile.
No saved logarithm/epsilon table, old global config, or symmetric proof
archive is read.  Optional --compare-* flags are one-time migration audits.
"""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
from time import monotonic
import argparse
import json
import os
import tempfile
from flint import arb, fmpq, ctx

BITS, LOG_N, EPS_N = 48, 16384, 1024
SCALE = 1 << BITS
COMPONENTS = [
    ("cap", "1", None, "166466453839/1000000000000"),
    ("old", "37037037037/62500000000", "3/4", "3055240759/62500000000"),
    ("old", "125/256", "4/5", "44226298717/1000000000000"),
    ("old", "258117479171/1000000000000", "9/10", "8208622709/250000000000"),
    ("old", "196045038119/1000000000000", "37/40", "21914693827/1000000000000"),
    ("old", "1358461287/500000000000", "999/1000", "254146426697/1000000000000"),
    ("cap", "1/2000", None, "15666000259/100000000000"),
    ("cap", "1/2", None, "5497355627/20000000000"),
]
PROXIES = [
    ("0", "0", "395494177/1000000000", "10002401/10000000"),
    ("4", "0", "2/5", "100032357/100000000"),
    ("4", "0", "1/2", "25657963/25000000"),
    ("4", "1/10", "9/10", "119092969/100000000"),
    ("4", "1/10", "1", "4949053/4000000"),
]
PLANAR = [("0", "27140539/100000000", "36074759/50000000")]


def aq(q):
    q = Q(q)
    return arb(fmpq(q.numerator, q.denominator))


def endpoint(value, upper=False):
    point = value.upper() if upper else value.lower()
    mantissa, exponent = point.man_exp()
    return Q(int(mantissa)) * Q(2)**int(exponent)


def rational_bounds(q):
    q = Q(q)*SCALE
    return q.numerator//q.denominator, -((-q.numerator)//q.denominator)


def bounds(value):
    lo = rational_bounds(endpoint(value))[0]
    hi = rational_bounds(endpoint(value, True))[1]
    assert aq(Q(lo,SCALE)) <= value <= aq(Q(hi,SCALE))
    return lo, hi


def epsilon(r):
    r = aq(r)
    ell = (1/r).log1p()
    return r+((1/(r*ell)).log()-1)/ell


def atomic(path, data):
    assert not path.exists(), f"Refusing to overwrite {path}"
    with tempfile.NamedTemporaryFile("wb", dir=path.parent, delete=False) as f:
        tmp = Path(f.name)
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp,path)
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def parse_tables(path):
    # Used only by the optional audit; every generated runtime value is
    # computed independently before this old text is opened.
    toks=iter(path.read_text().split())
    def tag(s):
        assert next(toks)==s
    def integer(): return int(next(toks))
    def rows(n): return [(integer(),integer()) for _ in range(n)]
    tag('P'); assert integer()==BITS
    tag('LOG_N'); assert integer()==LOG_N
    tag('LOG'); logs=rows(integer())
    tag('EPS_N'); assert integer()==EPS_N
    tag('EPS_BLOCKS'); n=integer()
    blocks=[]
    for e in range(-10,10):
        tag('EPS'); assert integer()==e
        blocks.append(rows(integer()))
    assert n==20
    tag('CONSTANTS'); n=integer()
    constants={}
    for _ in range(n):
        name=next(toks);constants[name]=(integer(),integer())
    tag('END')
    assert next(toks,None) is None
    return logs,blocks,constants


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--private-profile',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--compare-input',type=Path)
    ap.add_argument('--compare-tables',type=Path)
    args=ap.parse_args()
    started=monotonic();ctx.prec=192
    args.output_dir.mkdir(parents=True,exist_ok=True)
    private_raw=args.private_profile.read_bytes()
    profile=json.loads(private_raw)
    canonical=[(x['kind'],Q(x['theta']),Q(x['A']) if 'A' in x else None,Q(x['weight']))
               for x in profile['components']]
    expected=[(kind,Q(theta),Q(a) if a else None,Q(w)) for kind,theta,a,w in COMPONENTS]
    assert canonical==expected
    assert sum(x[-1] for x in expected)==1
    rows=profile['rows']
    assert len(rows)==4001
    assert all(Q(x['m'])==Q(j,1000) for j,x in enumerate(rows))
    values=[Q(x['lower_G']) for x in rows]
    assert all(0<a<=b for a,b in zip(values,values[1:]))
    lines=['P 48','PROXIES 5']
    for row in PROXIES:
        lines.append(' '.join(str(v) for q in row for v in rational_bounds(q)))
    lines.append('PLANAR 1')
    for row in PLANAR:
        lines.append(' '.join(str(v) for q in row for v in rational_bounds(q)))
    lines.append('PROFILE 4001')
    lines.extend(f'{lo} {hi}' for lo,hi in map(rational_bounds,values))
    lines.append('END')
    input_raw=('\n'.join(lines)+'\n').encode()

    z=arb(-1).exp();beta=1-z;ell=arb(2).log()
    rstar=(1/ell).log()/ell;epsstar=1+rstar-1/ell
    named={'Z':z,'beta':beta,'kappa':1-2*z,'L':ell,'rstar':rstar,
           'zeta':z-epsstar,'home_join':beta+rstar,
           'b':aq(Q(9,16)),'k_M':aq(Q(5,16)),'eta':aq(Q(3,5)),
           'c0':aq(Q(3999999,10000000))}
    constants={name:bounds(value) for name,value in named.items()}
    logs=[bounds(aq(Q(LOG_N+i,LOG_N)).log()) for i in range(LOG_N+1)]
    assert logs[0]==(0,0) and logs[-1]==constants['L']
    assert all(a[1]<b[0] for a,b in zip(logs,logs[1:]))
    blocks=[]
    for e in range(-10,10):
        block=[bounds(epsilon(Q(2)**e*Q(EPS_N+i,EPS_N))) for i in range(EPS_N+1)]
        assert all(0<lo<=hi<SCALE for lo,hi in block)
        assert all(a[0]>b[1] for a,b in zip(block,block[1:]))
        if blocks: assert blocks[-1][-1]==block[0]
        blocks.append(block)
    tablelines=['P 48',f'LOG_N {LOG_N}',f'LOG {len(logs)}']
    tablelines.extend(f'{lo} {hi}' for lo,hi in logs)
    tablelines.extend([f'EPS_N {EPS_N}','EPS_BLOCKS 20'])
    for e,block in zip(range(-10,10),blocks):
        tablelines.append(f'EPS {e} {len(block)}')
        tablelines.extend(f'{lo} {hi}' for lo,hi in block)
    tablelines.append(f'CONSTANTS {len(constants)}')
    tablelines.extend(f'{name} {lo} {hi}' for name,(lo,hi) in constants.items())
    tablelines.append('END')
    table_raw=('\n'.join(tablelines)+'\n').encode()
    comparisons={}
    if args.compare_input:
        old=args.compare_input.read_bytes()
        assert old==input_raw
        comparisons['rational_input']={'status':'BYTE_IDENTICAL','sha256':sha256(old).hexdigest()}
    if args.compare_tables:
        oldlogs,oldblocks,oldconst=parse_tables(args.compare_tables)
        assert logs==oldlogs and blocks==oldblocks
        assert all(oldconst[name]==pair for name,pair in constants.items())
        comparisons['used_table_entries']={'status':'ALL_IDENTICAL','nodes':len(logs)+sum(map(len,blocks)),
                                           'constants':len(constants),
                                           'old_file_sha256':sha256(args.compare_tables.read_bytes()).hexdigest()}
    report={'status':'PASS_LEAN_RUNTIME_DATA','precision_bits':192,'denominator_power':BITS,
            'log_nodes':len(logs),'epsilon_nodes':sum(map(len,blocks)),'constants':len(constants),
            'maximum_endpoint_width':max(hi-lo for lo,hi in logs+[x for b in blocks for x in b]+list(constants.values())),
            'private_nodes':4001,'proxies':5,'planar':1,'runtime_score_terms':8,
            'private_profile_sha256':sha256(private_raw).hexdigest(),
            'generator_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'table_sha256':sha256(table_raw).hexdigest(),'input_sha256':sha256(input_raw).hexdigest(),
            'comparisons':comparisons,'seconds':monotonic()-started,
            'scope':'Directed data generation, not a replay of the local score/profile certificates or global cover.'}
    atomic(args.output_dir/'tables48.txt',table_raw)
    atomic(args.output_dir/'input48.txt',input_raw)
    atomic(args.output_dir/'runtime_data_report.json',(json.dumps(report,indent=2)+'\n').encode())
    print(report['status'],report['log_nodes']+report['epsilon_nodes'],'transcendental nodes;',len(constants),'constants')
    print('input',report['input_sha256'])
    print('table',report['table_sha256'])
    print('seconds',report['seconds'])


if __name__=='__main__':
    main()
