#!/usr/bin/env python3
"""Recompute every numerical obligation of the accompanying simplified proof.

Requires Python 3, python-flint, and g++ (C++17).  No saved PASS, old package,
cached table, optimizer maximum, or previous cover is accepted as evidence.
The mathematical reductions and oracle algorithm are proved in paper/main.tex;
this driver is a complete computational replay, not a proof assistant.
"""
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
from time import monotonic
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import threading

ROOT = Path(__file__).resolve().parent
C = ROOT / 'certificate'
SCALE = 1 << 48
U_LOWER = (Q('0.6905799') * SCALE).__floor__()
YMAX = {1: '128', 2: '32', 3: '6', 4: '2'}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def atomic(path, obj):
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as f:
        temp = Path(f.name)
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def run_logged(command, log, env):
    start = monotonic()
    with log.open('wb') as f:
        child = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, env=env)
        f.flush()
        os.fsync(f.fileno())
    assert child.returncode == 0, (command, child.returncode, log.read_text()[-4000:])
    return {'command': command, 'exit_code': 0, 'seconds': monotonic()-start,
            'log_sha256': digest(log)}


def terminal(log, h):
    lines = log.read_text().splitlines()
    assert not any('INCOMPLETE' in s or 'Assertion' in s for s in lines)
    ends = [s for s in lines if s.startswith('COMPLETE ')]
    pattern = (r'COMPLETE B_SIMPLE h (\d+) boxes (\d+) accepted (\d+) excluded (\d+) '
               r'analytic_band (\d+) roots (\d+) depth (\d+) highest_upper (-?\d+) / (\d+) '
               r'threshold_lower (\d+) / (\d+) seconds ([0-9.eE+-]+)')
    assert len(ends) == 1
    m = re.fullmatch(pattern, ends[0])
    assert m, ends
    H, boxes, accepted, excluded, band, roots, depth, upper, scale, threshold, scale2 = map(int, m.groups()[:-1])
    assert H == h and roots == 2 and scale == scale2 == SCALE
    assert boxes == 2*(accepted+excluded)-roots
    assert 0 <= band < accepted and 0 <= depth < 180
    assert upper < threshold == U_LOWER
    starts = [s for s in lines if s.startswith('START ')]
    assert len(starts) == 1 and starts[0].startswith(f'START B_SIMPLE h {h} ')
    st = starts[0].split()
    assert int(st[st.index('scale')+1]) == SCALE
    assert int(st[st.index('roots')+1]) == 2
    assert int(st[st.index('U')+1]) == U_LOWER
    da = Q(21, 20)*Q('1.00032357')**2*SCALE
    dal, dah = map(int, st[st.index('A_D')+1:st.index('A_D')+3])
    assert dal <= da <= dah
    initial = [s for s in lines if s.startswith('BOX depth 0 ')]
    assert len(initial) == 1
    tokens = initial[0].split()
    tmax = (Q('2.500000626')*SCALE).__ceil__()
    assert list(map(int, tokens[4:6])) == [0, tmax]
    assert list(map(int, tokens[7:9])) == [0, 21*SCALE]
    for j in range(h):
        idx = tokens.index(f'y{j}')
        assert list(map(int, tokens[idx+1:idx+3])) == [0, int(YMAX[h])*SCALE]
        idx = tokens.index(f'S{j}')
        assert list(map(int, tokens[idx+1:idx+3])) == [0, tmax]
    return {'h': h, 'boxes': boxes, 'accepted': accepted, 'excluded': excluded,
            'analytic_band': band, 'roots': roots, 'depth': depth, 'highest_upper': upper,
            'denominator': SCALE, 'threshold_lower': threshold,
            'directed_numerical_gap': str(Q(threshold-upper, SCALE)),
            'analytic_gap': '1/1000000', 'terminal': ends[0]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path, required=True, help='Must not already exist')
    ap.add_argument('--jobs', type=int, default=4, choices=range(1, 5))
    args = ap.parse_args()
    assert sys.flags.optimize == 0, 'Run Python without -O'
    dest = args.output_dir.resolve()
    assert not dest.exists(), f'Fresh directory required: {dest}'
    dest.mkdir(parents=True)
    logs = dest/'logs'
    logs.mkdir()
    env = dict(os.environ)
    env.pop('PYTHONOPTIMIZE', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    sources = [Path(__file__).resolve()] + sorted((ROOT/'paper').glob('*.tex'))
    sources += sorted(p for p in C.rglob('*') if p.suffix in ('.py', '.cpp', '.json'))
    original = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    report = {'status': 'RUNNING', 'fresh_replay': True, 'started_utc': utc(),
              'source_sha256': original, 'gates': [], 'bounded_runs': []}
    start = monotonic()
    report_path = dest/'replay_report.json'
    atomic(report_path, report)
    profile, parameters = C/'local/private_profile.json', C/'local/parameters.json'

    def gate(name, command, output, expected):
        assert not output.exists()
        result = run_logged(command, logs/(name+'.log'), env)
        data = json.loads(output.read_text())
        assert data['status'] == expected, (name, data.get('status'))
        result.update({'name': name, 'status': expected, 'report': str(output.relative_to(dest)),
                       'report_sha256': digest(output), 'log': f'logs/{name}.log'})
        report['gates'].append(result)
        print('GATE_COMPLETE', name, round(result['seconds'], 3), 'seconds', flush=True)
        atomic(report_path, report)
        return data

    def py(name):
        return [sys.executable, str(C/name)]

    try:
        local = gate('local', py('local/replay_local.py') + ['--parameters', str(parameters),
            '--private-profile', str(profile), '--output', str(dest/'local.json'),
            '--runtime-generator', str(C/'generate_runtime_data.py')], dest/'local.json',
            'PASS_MINIMAL_LOCAL_CERTIFICATES')
        assert local['runtime_handoff_checked'] and local['private_profile']['points'] == 4001
        assert local['unconditional']['checks'] == 727 and local['conditional']['checks'] == 37304
        assert local['joint_planar']['vertices'] == 64
        minor = gate('minorants', py('check_minorants.py') + [str(parameters), str(dest/'minorants.json')],
                     dest/'minorants.json', 'PASS_COMPLETION_AND_MINORANTS')
        assert minor['checks'] == 51 and minor['high_cover']['intervals'] == 467
        home = gate('home', py('check_home.py') + [str(profile), str(dest/'home.json')],
                    dest/'home.json', 'PASS_COMPACT_GLOBAL_HOME')
        assert home['profile_intervals'] == 4000
        gate('A', py('check_scalar_A.py') + ['--output', str(dest/'A.json')],
             dest/'A.json', 'PASS_CONSOLIDATED_SCALAR_A')
        band = gate('band', py('check_scalar_band.py') + ['--output', str(dest/'band.json')],
                    dest/'band.json', 'PASS_CONSOLIDATED_SCALAR_BAND')
        assert band['check_count'] == 80 and band['counts'] == [1, 2, 3, 4]
        assert band['positive_count_dependencies'] == [] and Q(band['uniform_reserve']) == Q(1, 10**6)
        assert band['domain'] == {'t': ['1199/10000', '3001/10000'], 'D': ['0', '23/20']}
        gate('finite', py('check_finite.py') + ['--output', str(dest/'finite.json')],
             dest/'finite.json', 'PASS_GRAPH_AND_FINITE_INTERFACE')
        gate('lower_bound', py('check_lower_bound.py') + ['--output', str(dest/'lower_bound.json')],
             dest/'lower_bound.json', 'PASS_SIMPLIFIED_RATIONAL_FRAMEWORK_LOWER_BOUND')
        tail = gate('tail', py('check_tail.py') + [str(profile), str(dest/'tail.json')],
                    dest/'tail.json', 'PASS_TAIL_AND_ALL_COUNTS')
        assert Q(tail['strict_gap']) > 0
        assert {int(h): str(v['cutoff']) for h, v in tail['deletions'].items()} == {h: YMAX[h] for h in (2, 3, 4)}
        gate('menu', py('check_menu_reductions.py') + ['--components', str(profile),
            '--source', str(C/'verify_B_simple.cpp'), '--output', str(dest/'menu.json')],
            dest/'menu.json', 'PASS_EXACT_MENU_COMPRESSION_AND_B_REDUNDANCY')
        data = gate('runtime', py('generate_runtime_data.py') + ['--private-profile', str(profile),
            '--output-dir', str(dest/'data')], dest/'data/runtime_data_report.json', 'PASS_LEAN_RUNTIME_DATA')
        assert data['runtime_score_terms'] == 8 and data['private_nodes'] == 4001
        table, inp, executable = dest/'data/tables48.txt', dest/'data/input48.txt', dest/'verify_B'
        assert digest(table) == data['table_sha256'] and digest(inp) == data['input_sha256']
        source = C/'verify_B_simple.cpp'
        assert 'if(b.t.l>=rat(1199,10000).h && b.t.h<=rat(3001,10000).l && b.D.h<=rat(23,20).l)' in source.read_text()
        command = ['g++', '-std=c++17', '-O3', '-Wall', '-Wextra', str(source), '-o', str(executable)]
        build = run_logged(command, logs/'compile.log', env)
        build.update({'source_sha256': digest(source), 'executable_sha256': digest(executable),
            'compiler': subprocess.run(['g++', '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]})
        atomic(dest/'build_receipt.json', build)
        gate('primitives', py('audit_simple_B.py') + ['--source', str(source), '--executable', str(executable),
            '--tables', str(table), '--input', str(inp), '--output', str(dest/'primitives.json')],
            dest/'primitives.json', 'PASS_SIMPLE_B_INDEPENDENT_PRIMITIVES')
        fixed = {p: digest(p) for p in (source, table, inp, executable)}
        receipts, lock = [], threading.Lock()

        def launch(h):
            before = {str(p): digest(p) for p in fixed}
            assert before == {str(p): v for p, v in fixed.items()}
            command = [str(executable), str(table), str(inp), str(h), '30000000',
                       '0', '2.500000626', '0', '21', YMAX[h]]
            receipt = {'h': h, 'started_utc': utc(), 'status': 'RUNNING',
                       'freshly_executed': True, 'command': command, 'prelaunch_sha256': before}
            with lock:
                receipts.append(receipt)
                atomic(dest/'launch_receipts.json', receipts)
            log = logs/f'h{h}.log'
            run = run_logged(command, log, env)
            after = {str(p): digest(p) for p in fixed}
            assert before == after
            result = terminal(log, h)
            result.update(run)
            result.update({'freshly_executed': True, 'log': f'logs/h{h}.log'})
            with lock:
                receipt.update({'status': 'COMPLETE', 'exit_code': 0, 'finished_utc': utc(),
                    'postrun_sha256': after, 'log_sha256': digest(log)})
                atomic(dest/'launch_receipts.json', receipts)
            return result

        print('FRESH_COVERS_STARTED', dest, flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending = {pool.submit(launch, h): h for h in range(1, 5)}
            while pending:
                done, _ = wait(pending, timeout=20, return_when=FIRST_COMPLETED)
                for f in done:
                    h = pending.pop(f)
                    result = f.result()
                    report['bounded_runs'].append(result)
                    print('COVER_COMPLETE', h, result['boxes'], 'boxes', round(result['seconds'], 3), 'seconds', flush=True)
                if not done:
                    print('REPLAY_RUNNING', round(monotonic()-start, 1), 'seconds; pending counts', sorted(pending.values()), flush=True)
                atomic(report_path, report)
        assert sorted(x['h'] for x in report['bounded_runs']) == [1, 2, 3, 4]
        assert original == {str(p.relative_to(ROOT)): digest(p) for p in sources}
        report.update({'status': 'PASS_COMPLETE_SIMPLIFIED_TIGHT_PROOF_REPLAY',
            'finished_utc': utc(), 'seconds': monotonic()-start,
            'bounded_runs': sorted(report['bounded_runs'], key=lambda x: x['h']),
            'total_boxes': sum(x['boxes'] for x in report['bounded_runs']),
            'all_proof_and_certificate_inputs_unchanged': True,
            'build_receipt_sha256': digest(dest/'build_receipt.json'),
            'launch_receipts_sha256': digest(dest/'launch_receipts.json'),
            'scope': 'All finite computational obligations of the accompanying mathematical proof were freshly executed. General lemmas, rounding and oracle complexity are proved in the manuscript; this is not formal proof-assistant verification.'})
        atomic(report_path, report)
        print(report['status'], report['total_boxes'], 'boxes;', round(report['seconds'], 3), 'seconds', flush=True)
    except BaseException as exc:
        report.update({'status': 'FAIL_OR_INCOMPLETE', 'error': repr(exc), 'seconds': monotonic()-start})
        atomic(report_path, report)
        raise


if __name__ == '__main__':
    main()
