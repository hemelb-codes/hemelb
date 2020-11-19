# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import io
import pytest
import numpy as np

import parareal
from parareal.interval import Interval
from parareal.interfaces import hm

@pytest.fixture
def clean_luigi(monkeypatch):
    '''Luigi has singletons for various configuration objects and an
    instance cache of horror. Here we monkey patch those back to the
    initial state so it's like we have just done `import luigi`
    '''
    import luigi

    monkeypatch.setattr(luigi.configuration.LuigiConfigParser,
                            '_instance',
                            None)
    monkeypatch.setattr(luigi.task_register.Register,
                            '_Register__instance_cache',
                            {})
    monkeypatch.setattr(luigi.cmdline_parser.CmdlineParser,
                            '_instance',
                            None)

    monkeypatch.setattr(luigi.interface.setup_interface_logging,
                            'has_run',
                            False,
                            raising=False)
    return luigi

def set_params(clean_luigi, params):
    cfg = clean_luigi.configuration.LuigiConfigParser.instance()
    
    cfg.set('core', 'workers', '1')
    cfg.set('core', 'no_configure_logging', 'True')
    cfg.set('core', 'local_scheduler', 'True')
    
    params.to_file('tmp.yml')
    cfg.set('parareal', 'filename', 'tmp.yml')
    return cfg

def test_simple_runs(tmpdir, clean_luigi):
    params = parareal.PararealParameters(
        {'b': 0.1, 'omega': 1.0},
        {'x': 1.0, 'v': 0.0},
        Interval(2.0 * np.pi/1000,0.0, 1000),
        2, 4)
    
    with tmpdir.as_cwd() as old:
        set_params(clean_luigi, params)
        tsk = params.final_task()
        assert not tsk.output().exists()
        ran_ok = clean_luigi.build([tsk])
        assert ran_ok
        assert tsk.output().exists()

def merge_trajectories(traj_lst):
    chunk_lst = []
    n_rows = 0
    for t in traj_lst:
        chunk = np.loadtxt(t.path)
        try:
            last = chunk_lst[-1]
            assert np.all(last[-1] == chunk[0])
            chunk = chunk[1:]
        except IndexError:
            # we must be first
            pass
        chunk_lst.append(chunk)
        continue
    ans = np.concatenate(chunk_lst, axis=0)
    np.savetxt('parareal.txt', ans)
    return ans

def run_serial(parareal_params):
    prob = hm.Problem(parareal_params.sim_params,
                          parareal_params.initial_conditions,
                          parareal_params.time,
                          'fine')
    fn = 'serial.txt'
    prob.run(fn)
    return np.loadtxt(fn)
    
    
def test_44_equiv_serial(tmpdir, clean_luigi, plot_mode=False):
    i_max = 4
    
    params = parareal.PararealParameters(
        {'b': 0.1, 'omega': 1.0},
        {'x': 1.0, 'v': 0.0},
        Interval(2.0 * np.pi/1000,0.0, 1000),
        i_max, i_max)
    
    with tmpdir.as_cwd() as old:
        set_params(clean_luigi, params)
        tsk = params.final_task()
        
        ran_ok = clean_luigi.build([tsk])
        assert ran_ok
        
        assert tsk.output().exists()

        overall = merge_trajectories(
            parareal.Trajectory(parareal.Resolution.Fine, i, i)
            for i in range(i_max)
            )
        serial = run_serial(params)

        delta = serial - overall
        
        dt, dx, dv = delta.transpose()
        
        if not plot_mode:
            # We have to choose these tolerances so skip this to get
            # the data out and plotted.
            assert np.all(dt == 0)

            e2 = np.sqrt(dx**2 + dv**2)
            assert np.all(np.absolute(e2) < 5e-4)
            
