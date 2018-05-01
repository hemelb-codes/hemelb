import os.path
import enum
from functools import wraps

import yaml
import numpy as np
import luigi

from interval import Interval
import hm

import pdb

class Resolution(enum.Enum):
    """Flag which resolution we're at"""
    Coarse = 0
    Fine = 1
    
    @property
    def other(self):
        if self == Resolution.Coarse:
            return Resolution.Fine
        else:
            return Resolution.Coarse
        
    pass


class PararealParameters(object):
    """This class is responsible for describing the overall set of
    simulations and performing some set up on behalf of the actual luigi
    tasks.
    """

    def __init__(self, params, initial_conditions, time, num_parareal_iters, num_time_slices):
        """
        params - parameters needed to do a time integration.
        
        initial_conditions - overall ICs

        time - overall time interval to be integrated

        num_parareal_iters - number of PR iters
        
        num_time_slices - how many sub intervals to split the overall interval into
        
        """

        self.params = params
        self.initial_conditions = initial_conditions
        self.time = time
        self.num_parareal_iters = num_parareal_iters
        self.num_time_slices = num_time_slices
        
        self.sub_times = time.subdivide(num_time_slices)
        
        return

    def write_icond_input(self, input_target):
        self.write_input(input_target, self.initial_conditions)
        
    def write_input(self, input_target, ic):
        sub_ival = self.sub_times[input_target.i]
        res = {
            Resolution.Fine: 'fine',
            Resolution.Coarse: 'coarse'
            }[input_target.res]
            
        state = {
            'params': self.params,
            'initial': ic,
            'time': sub_ival.to_dict(),
            'solver': input_target.res.name.lower()
            }
        with input_target.open('w') as f:
            yaml.dump(state, stream=f)

    def final_task(self):
        return ParaUpdate(res=Resolution.Fine, i=self.num_time_slices, k=self.num_parareal_iters)
    
    pass

# luigi might run tasks in a new process, make sure the PP are available
# to it. Yes this is horrible - should probably be a task and target to
# load from a config file.
_master = None
def get_master():
    global _master
    
    if _master is None:
        num_time_slices = 4
        num_parareal_iters = 2

        params = {'b': 0.1, 'omega': 1.0}
        initial = {
                'x': 1.0,
                'v': 0.0
            }
        time = Interval(2.0 * np.pi/1000, 0, 1000)

        _master = PararealParameters(params, initial, time, num_parareal_iters, num_time_slices)
    return _master


class ParaRealTarget(luigi.LocalTarget):
    """A local file pattern, specialised based on resolution, i, k.

    Subclass should set _pattern before this constructor is called.
    """
    def __hash__(self):
        return hash(self.path)
    def __eq__(self, other):
        return self.path == other.path

    def __init__(self, res, i, k):
        self.res = res
        if i == 0:
            k = 0
            pass
        
        self.i = i
        self.k = k
        resolution = res.name.lower()
        path = self._pattern.format(**locals())
        super(ParaRealTarget, self).__init__(path)
        return

    pass

def check_producer(prod_func):
    @wraps(prod_func)
    def wrapper(cls, res, i, k):
        tsk = prod_func(cls, res, i, k)
        if tsk.output() != cls(res, i, k):
            pdb.set_trace()
        return tsk
    return wrapper

class Input(ParaRealTarget):
    """Represent an input file, ready to be run.
    """
    _pattern ='{resolution}_{i:04d}_{k:02d}/problem.yml'

    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        assert i >= 0
        assert k >= 0
        assert k <= i
        
        if i == 0:
            return InitialConditionMaker(res)
        else:
            return InputMaker(res=res, i=i, k=k)
    pass

class Trajectory(ParaRealTarget):
    _pattern = '{resolution}_{i:04d}_{k:02d}/results.txt'
    
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        return Solver(res=res, i=i, k=k)

    def get_last(self):
        data = np.loadtxt(self.open())
        return data[-1].tolist()
    
    pass

class State(ParaRealTarget):
    """A state is a yaml file containing t, x, v (floats)
    """

    def load(self):
        with self.open() as reader:
            return yaml.load(reader)

    def save(self, obj):
        with self.open('w') as writer:
            yaml.dump(obj, stream=writer)


    pass

class y(State):
    _pattern = 'state/y_{resolution}_{i:04d}_{k:02d}'
            
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        assert i >= 0
        assert k >= 0
        assert k <= i

        ans = None
        if i == 0:
            raise ValueError('Should probably never request y[i=0] explicitly')
        else:
            if k == 0:
                # y[i, k=0] = G(y[i-1, k=0])
                return Linker(res=res, i=i, k=k)
            
            elif k == i:
                # y[i,k=i] = F(y(i=i-1, k=i-1))
                return Linker(res=res, i=i, k=k)
            
            else:
                return ParaUpdate(res=res, i=i, k=k)
    
    # def __str__(self):
    #     return "y[i={}, k={}]".format(self.i, self.k)
    pass

class DeltaGy(State):
    _pattern = 'state/DGy_{resolution}_{i:04d}_{k:02d}'
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src_cls = {
            Resolution.Fine: DiffRefinement,
            Resolution.Coarse: CoarseDiffer
            }[res]
        return src_cls( i=i, k=k)
    pass

class Gy(State):
    _pattern = 'state/Gy_{resolution}_{i:04d}_{k:02d}'
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src_cls = {
            Resolution.Fine: SolutionRefinement,
            Resolution.Coarse: SolutionExtractor
            }[res]
        return src_cls(res=Resolution.Coarse, i=i, k=k)
        
    # def __str__(self):
    #     return "G(y[i={}, k={}])".format(self.i, self.k)
    pass

class Fy(State):
    _pattern = 'state/Fy_{resolution}_{i:04d}_{k:02d}'
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src = SolutionExtractor(res=Resolution.Fine, i=i, k=k)
        
        if res == Resolution.Coarse:
            return Coarsening(src)
        
        return src
    
    # def __str__(self):
    #     return "F(y[i={}, k={}])".format(self.i, self.k)

    pass

class InitialConditionMaker(luigi.ExternalTask):
    # Which resolution 
    res = luigi.EnumParameter(enum=Resolution)
    
    def output(self):
        return Input(self.res, 0, 0)
    
    def run(self):
        mst = get_master()
        mst.write_icond_input(self.output())
        return
    pass


class Op(luigi.Task):    
    # time slice ID
    i = luigi.IntParameter()
    # parareal ID
    k = luigi.IntParameter()
    # Which resolution 
    res = luigi.EnumParameter(enum=Resolution)
    
    pass

class InputMaker(Op):
    def output(self):
        return Input(self.res, self.i, self.k)
    def requires(self):
        return y.producer(self.res, self.i, self.k)
    def run(self):
        state = self.input().load()
        ic = {'x': state['x'], 'v': state['v']}
        mst = get_master()
        mst.write_input(self.output(), ic)
        return
    pass

class Solver(Op):
    """Do one iteration of the basic solvers.
    
    Advances time from i -> i+1
    """
    
    def output(self):
        return Trajectory(self.res, self.i, self.k)
    
    def requires(self):
        return Input.producer(self.res, self.i, self.k)
    
    def run(self):
        prob = hm.Problem.from_file(self.input().path)
        prob.run(self.output().path)
        return
    
    pass

class SolutionExtractor(Op):
    def output(self):
        cls = {Resolution.Coarse: Gy, Resolution.Fine: Fy}[self.res]
        return cls(self.res, self.i, self.k)
    def requires(self):
        return Trajectory.producer(self.res, self.i, self.k)
    def run(self):
        t, x, v = self.input().get_last()
        self.output().save({'t': t, 'x': x, 'v': v})
        return
    pass

class Linker(Op):
    def requires(self):
        if self.k == 0:
            # y[i, k=0] = G(y[i-1, k=0])
            if self.res == Resolution.Coarse:
                return SolutionExtractor(res=Resolution.Coarse, i=i-1, k=0)
            else:
                return SolutionRefinement(src)
        
        elif self.k == self.i:
            # y[i,k=i] = F(y(i=i-1, k=i-1))
            if self.res == Resolution.Fine:
                return SolutionExtractor(res=Resolution.Fine, i=i-1, k=i-1)
            else:
                return Coarsening(src)
        else:
            raise ValueError("should never get here")
        
    def output(self):
        return y(self.res, self.i, self.k)
    
    def run(self):
        os.link(self.input().path, self.output().path)
        return

class Coarsening(luigi.Task):
    # TODO: make this work like Refinement and subclasses
    source = luigi.TaskParameter()

    def requires(self):
        return self.source
    
    def output(self):
        instate = self.input()
        assert instate.res == Resolution.Fine
        return type(instate)(Resolution.Coarse, instate.i, instate.k)

    def run(self):
        self.output().save(self.input().load())
        return
    pass

class Refinement(luigi.Task):
    i = luigi.IntParameter()
    k = luigi.IntParameter()
    def requires(self):
        return self.refinee_cls.producer(Resolution.Coarse, self.i, self.k)
    def output(self):
        return self.refinee_cls(Resolution.Fine, self.i, self.k)

    def run(self):
        self.output().save(self.input().load())
        return
    
class SolutionRefinement(Refinement):
    refinee_cls = Gy    
    pass

class DiffRefinement(Refinement):
    refinee_cls = DeltaGy
    pass

class CoarseDiffer(luigi.Task):
    i = luigi.IntParameter()
    k = luigi.IntParameter()
    
    def requires(self):
        return [
            Gy.producer(Resolution.Coarse, self.i, self.k),
            Gy.producer(Resolution.Coarse, self.i, self.k-1)
            ]
    
    def output(self):
        return DeltaGy(Resolution.Coarse, self.i, self.k)

    def run(self):
        gk, gk_1 = [i.load() for i in self.input()]
        # TODO: check time is equal to times[self.i + 1]
        assert gk['t'] == gk_1['t']

        state = {
            't': gk['t'],
            'x': gk['x'] - gk_1['x'],
            'v': gk['v'] - gk_1['v'],
            }
        self.output().save(state)
        return
    pass

class ParaUpdate(Op):
    def output(self):
        return y(self.res, self.i, self.k)
    
    def requires(self):
        return [
            DeltaGy.producer(self.res, self.i-1, self.k),
            Fy.producer(self.res, self.i-1, self.k-1)
            ]
    def run(self):
        DGy, Fy = [i.load() for i in self.input()]
        # TODO: check time is equal to times[self.i + 1]
        assert DGy['t'] == Fy['t']

        state = {
            't': Fy['t'],
            'x': Fy['x'] + DGy['x'],
            'v': Fy['v'] + DGy['v'],
            }
        self.output().save(state)
        return
    pass



if __name__ == '__main__':
    tsk = y.producer(Resolution.Fine, 4, 2)
    # #tsk = Solver(res=Resolution.Coarse, i=3, k=0)
    luigi.build([tsk], local_scheduler=True)
