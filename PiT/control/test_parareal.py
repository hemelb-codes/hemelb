import pytest
import parareal
import luigi

Coarse = parareal.Resolution.Coarse
Fine = parareal.Resolution.Fine

@pytest.mark.parametrize("i", [1,2,5])
def test_y0_i_is_made_by_coarse(i):
    t0 = parareal.y.producer(Coarse, i, 0)
    assert isinstance(t0, parareal.Assignment)
    
    t1 = t0.requires()
    assert isinstance(t1, parareal.SolutionExtractor)
    
    t2 = t1.requires()
    assert isinstance(t2, parareal.Solver)
    assert t2.res == Coarse
    return

@pytest.mark.parametrize("i", [1,2,5])
def test_yi_i_is_made_by_fine(i):
    t0 = parareal.y.producer(Fine, i, i)
    assert isinstance(t0, parareal.Assignment)
    
    t1 = t0.requires()
    assert isinstance(t1, parareal.SolutionExtractor)
    
    t2 = t1.requires()
    assert isinstance(t2, parareal.Solver)
    assert t2.res == Fine
    
    return

def gather_deps(tsk):
    deps = tsk.requires()
    if deps is None:
        return {}

    if isinstance(deps, luigi.Task):
        ans = gather_deps(deps)
        ans.add(deps)
        return ans

    ans = set()
    for d in deps:
        ans.update(gather_deps(d))
    return ans

def test_32_deps(monkeypatch):
    tsk = parareal.y.producer(Fine, 3, 2)
    deps = gather_deps(tsk)
    assert len(deps) == 34
    
