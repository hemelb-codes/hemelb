# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import pytest
import parareal
import luigi

Coarse = parareal.Resolution.Coarse
Fine = parareal.Resolution.Fine

@pytest.mark.parametrize("i", range(1,10))
def test_y0_i_is_made_by_coarse(i):
    t0 = parareal.y.producer(Coarse, i, 0)
    assert isinstance(t0, parareal.Assignment)
    
    t1 = t0.requires()
    assert isinstance(t1, parareal.SolutionExtractor)
    
    t2 = t1.requires()
    assert isinstance(t2, parareal.Solver)
    assert t2.res == Coarse
    return

@pytest.mark.parametrize("i", range(1,10))
def test_yi_i_is_made_by_fine(i):
    t0 = parareal.y.producer(Fine, i, i)
    assert isinstance(t0, parareal.Assignment)
    
    t1 = t0.requires()
    assert isinstance(t1, parareal.SolutionExtractor)
    
    t2 = t1.requires()
    assert isinstance(t2, parareal.Solver)
    assert t2.res == Fine
    
    return

def gather_deps(tsk, seen=None):
    if seen is None:
        seen = set()
    deps = tsk.requires()
    
    if not deps:
        return seen

    if isinstance(deps, luigi.Task):
        deps = [deps]
        
    for d in deps:
        if d in seen:
            continue
        seen.add(d)
        seen = gather_deps(d, seen)
        
    return seen

def test_gather_deps():
    tsk = luigi.Task()
    deps = gather_deps(tsk)
    assert deps == set()

    class T1(luigi.Task):
        def requires(self):
            return luigi.Task()
        pass
    class T2(luigi.Task):
        def requires(self):
            return T1()
        pass
    class T3(luigi.Task):
        def requires(self):
            return [T2(), T1()]
        pass

    tsk = T3()
    deps = gather_deps(tsk)
    assert deps == {T2(), T1(), luigi.Task()}
    
    
def test_graph_12():
    '''y^1_2 = G(y^1_1) + F(y^0_1) - G(y^0_1)
    
    y^1_1 = F(y^0_0)
    
    y^0_1 = G(y^0_0)

    But with all the conversion steps too!
    '''
        
    checked_tasks = set()
    def mark(t):
        assert t not in checked_tasks
        checked_tasks.add(t)
    
    tsk = parareal.y.producer(Coarse, 2, 1)
    mark(tsk)
    assert tsk == parareal.ParaUpdate(res=Coarse, i=2, k=1)

    # y^1_2 = G(y^1_1) + F(y^0_1) - G(y^0_1)
    dgy, fy = tsk.requires()
    mark(fy)
    assert fy == parareal.SolutionCoarsening(i=1,k=0)
    mark(dgy)
    assert dgy == parareal.CoarseDiffer(i=1, k=1)
    
    solex = fy.requires()
    mark(solex)
    assert solex == parareal.SolutionExtractor(res=Fine, i=1, k=0)
    
    solver = solex.requires()
    mark(solver)
    assert solver == parareal.Solver(res=Fine, i=1, k=0)
    
    im = solver.requires()
    mark(im)
    assert im == parareal.InputMaker(res=Fine, i=1, k=0)

    # y^0_1 = G(y^0_0)
    ass = im.requires()
    mark(ass)
    assert ass == parareal.Assignment(res=Fine, i=1, k=0)

    ref = ass.requires()
    mark(ref)
    assert ref == parareal.SolutionRefinement(i=0, k=0)

    solex = ref.requires()
    mark(solex)
    assert solex == parareal.SolutionExtractor(res=Coarse, i=0, k=0)

    solver = solex.requires()
    mark(solver)
    assert solver == parareal.Solver(res=Coarse, i=0, k=0)

    ic = solver.requires()
    mark(ic)
    assert ic == parareal.InitialConditionMaker(res=Coarse)

    # Finished this branch!
    assert not ic.requires()

    # Back to dgy (delta G(y))
    gy11, gy10 = dgy.requires()
    # G(y^1_1)
    mark(gy11)
    assert gy11 == parareal.SolutionExtractor(res=Coarse, i=1, k=1)
    # G(y^0_1)
    mark(gy10)
    assert gy10 == parareal.SolutionExtractor(res=Coarse, i=1, k=0)

    solver = gy11.requires()
    mark(solver)
    assert solver == parareal.Solver(res=Coarse, i=1, k=1)

    im = solver.requires()
    mark(im)
    assert im == parareal.InputMaker(res=Coarse, i=1, k=1)

    # y^1_1 = F(y^0_0)
    ass = im.requires()
    mark(ass)
    assert ass == parareal.Assignment(res=Coarse, i=1, k=1)

    coarsen = ass.requires()
    mark(coarsen)
    assert coarsen == parareal.SolutionCoarsening(i=0, k=0)

    solex = coarsen.requires()
    mark(solex)
    assert solex == parareal.SolutionExtractor(res=Fine, i=0, k=0)

    solver = solex.requires()
    mark(solver)
    assert solver == parareal.Solver(res=Fine, i=0, k=0)
    
    ic = solver.requires()
    mark(ic)
    assert ic == parareal.InitialConditionMaker(res=Fine)
    
    # Finished this branch!
    assert not ic.requires()
    
    # Back to G(y^0_1)
    solver = gy10.requires()
    mark(solver)
    assert solver == parareal.Solver(res=Coarse, i=1, k=0)

    im = solver.requires()
    mark(im)
    assert im == parareal.InputMaker(res=Coarse, i=1, k=0)

    # y^0_1 = G(y^0_0)
    # which we have seen above, but at res==Fine
    ass = im.requires()
    mark(ass)
    assert ass == parareal.Assignment(res=Coarse, i=1, k=0)

    # The solution extractor we have seen tho!
    solex = ass.requires()
    assert solex in checked_tasks
    
    # Thus ends the checking of the task graph!
    deps = gather_deps(tsk)
    deps.add(tsk)

    # check that we really do have exactly the right set of tasks with
    # no extra
    assert deps == checked_tasks

    # Check the outputs
    outs = {d.output() for d in deps}

    def check(op):
        if isinstance(op, set):
            assert op.issubset(outs)
            outs.difference_update(op)
        else:
            assert op in outs
            outs.remove(op)

    # Know final state is there
    check(parareal.y(Coarse, 2, 1))

    # We know which trajectories have to be computed
    Traj = parareal.Trajectory
    needed_trajs = {
        Traj(Coarse, 1,1), Traj(Fine, 1, 0), Traj(Coarse, 1, 0),
        Traj(Fine, 0, 0),
        Traj(Coarse, 0, 0)
        }
    check(needed_trajs)
        
    # We know that each trajectory has it's final state extracted
    def traj2state(t):
        cls = {Coarse: parareal.Gy, Fine: parareal.Fy}[t.res]
        return cls(t.res, t.i, t.k)
    
    check({
        traj2state(t) for t in needed_trajs
        })

    # We know that each solver needs an input
    check({
        parareal.Input(t.res, t.i, t.k)
        for t in needed_trajs
        })

    # Know that there is one paraupdate, so only one DeltaGy
    check(parareal.DeltaGy(Coarse, 1, 1))

    # Two coarsens
    check(parareal.Fy(Coarse, 1, 0))
    check(parareal.Fy(Coarse, 0, 0))
    # One refine
    check(parareal.Gy(Fine, 0, 0))
    # assign y^0_1 at both resolutions
    check(parareal.y(Fine, 1, 0))
    check(parareal.y(Coarse, 1, 0))
    # and y^1_1 at coarse
    check(parareal.y(Coarse, 1, 1))
    assert len(outs) == 0

@pytest.mark.parametrize("i", range(1,10))
def test_ii_equiv_serial(i):
    '''This should be just applying the fine solver repeatedly'''
    i = 4
    tsk = parareal.y.producer(Fine, i, i)

    while i:
        assert tsk.output() == parareal.y(Fine, i, i)
        assert tsk == parareal.Assignment(res=Fine, i=i, k=i)

        i -= 1
        tsk = tsk.requires()
        assert tsk == parareal.SolutionExtractor(res=Fine, i=i, k=i)

        tsk = tsk.requires()
        assert tsk == parareal.Solver(res=Fine, i=i, k=i)

        tsk = tsk.requires()
        if i == 0:
            assert tsk == parareal.InitialConditionMaker(res=Fine)
        else:
            assert tsk == parareal.InputMaker(res=Fine, i=i, k=i)    

        tsk = tsk.requires()
        
    return
