"""Peer equivalence of the two fibre engines, observed through the whole pipeline: `pcensus.py --engine pool` and
`--engine torch` must print the same JSON (all fields but wall time t). Fibres are sorted canonically before
materialization, so representative-dependent fields (transversal.lifted) are compared too.
The torch runs also set ETF_XCHECK=1: per parent and column, lift and children are compared with pcensus.lift/check."""
import json, os, subprocess, sys
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = [('K2', 3), ('K2', 5), ('pair_after', 3), ('pair_srcout', 3), ('pair_tgtin', 3)]


def run(q, p, *extra, env=None):
    out = subprocess.run([sys.executable, 'pcensus.py', q, str(p), *extra], cwd=HERE, capture_output=True, text=True,
                         env={**os.environ, **(env or {})}, check=True).stdout
    res = json.loads([l for l in out.splitlines() if l.startswith('{')][-1]); res.pop('t')
    return res


@pytest.mark.parametrize('q,p', CASES)
def test_pool_equals_torch(q, p):
    assert run(q, p, '--engine', 'pool', '--jobs', '2') == run(q, p, '--engine', 'torch', env={'ETF_XCHECK': '1'})


def test_checkpoint_crosses_engines(tmp_path):
    """A checkpoint written by one engine is resumed by the other (same format, same fibres)."""
    ck = str(tmp_path / 'ck.pkl')
    a = run('pair_after', 3, '--engine', 'pool', '--jobs', '2', '--checkpoint', ck)
    b = run('pair_after', 3, '--engine', 'torch', '--checkpoint', ck)          # all fibres resumed
    os.remove(ck)
    c = run('pair_after', 3, '--engine', 'torch', '--checkpoint', ck)
    d = run('pair_after', 3, '--engine', 'pool', '--jobs', '2', '--checkpoint', ck)
    assert a == b == c == d
