# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import os.path
import enum
from functools import wraps
from abc import abstractmethod
import logging

import yaml
import luigi

from .interval import Interval
from .iohelp import *

import sys
import subprocess
from contextlib import contextmanager
import os
import coarsen
import refine
import difference
import sum

numFineProcs = 4
numCoarseProcs = 3

coarseRankFileName = 'initialisation/coarseMpirank.xtr'
fineRankFileName = 'initialisation/fineMpirank.xtr'

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def run_hemelb(mpiTasks):
    sys.stdout.write('Running hemelb on ' + str(mpiTasks) + ' procs\n')
    runcommand = []
    runcommand.append("aprun")
    runcommand.append("-n")
    runcommand.append("{0}".format(mpiTasks))
    execName = "../hemelb"
    runcommand.append(execName)
    runcommand.append("-in")
    runcommand.append("config.xml")
    runcommand.append("-i")
    runcommand.append("1")
    runcommand.append("-ss")
    runcommand.append("1111")
    sys.stdout.write('Run command:')
    for item in runcommand:
        sys.stdout.write(' ' + item)
    sys.stdout.write('\n')
    p = subprocess.Popen(runcommand, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = p.communicate()
    sys.stdout.write(stdoutdata)
    sys.stderr.write(stderrdata)
    pass

class Resolution(enum.Enum):
    '''Flag which resolution we're at'''
    Coarse = 0
    Fine = 1
    
    @property
    def other(self):
        if self == Resolution.Coarse:
            return Resolution.Fine
        else:
            return Resolution.Coarse
        
    pass

class parareal(luigi.Config):
    '''Luigi might run tasks in a new process, so we have to make sure
    our PararealParameters instance is available to it. We make that
    file a parameter of this Config object, who's value is picked up by
    Luigi's standard mechanism.

    I suggest adding a section like

    [parareal]
    filename: foobar.yml
    
    to luigi.cfg in your run directory
    '''
    filename = luigi.Parameter()
    
    @property
    def data(self):
        try:
            # Ensure we don't use our custom getattr
            return super(parareal, self).__getattribute__('_data')
        except AttributeError:
            ans = self._data = PararealParameters.from_file(self.filename)
            return ans

    def __getattr__(self, attr):
        '''Delegate to the Parameters object'''
        return getattr(self.data, attr)
    
    pass


class PararealParameters(LoadableMixin, DumpableMixin):
    '''This class is responsible for describing the overall set of
    simulations and performing some set up on behalf of the actual luigi
    tasks.
    '''

    def __init__(self, sim_params, initial_conditions, time,
                     num_parareal_iters, num_time_slices):
        '''sim_params - parameters needed to do a time integration.
        
        initial_conditions - overall ICs

        time - overall time interval to be integrated

        num_parareal_iters - number of PR iters
        
        num_time_slices - how many sub intervals to split the overall
        interval into
        '''

        self.sim_params = sim_params
        self.initial_conditions = initial_conditions
        print('initial_conditions:')
        print(initial_conditions)
        self.time = time
        self.num_parareal_iters = num_parareal_iters
        self.num_time_slices = num_time_slices
        
        self.sub_times = time.subdivide(num_time_slices)
        
        return
    
    @classmethod
    def from_dict(cls, d):
        time = Interval.from_dict(d['time'])
        return cls(d['sim_params'], d['initial_conditions'], time,
                    d['num_parareal_iters'], d['num_time_slices'])

    def to_dict(self):
        return {
            'sim_params': self.sim_params,
            'initial_conditions': self.initial_conditions,
            'time': self.time.to_dict(),
            'num_parareal_iters': self.num_parareal_iters,
            'num_time_slices': self.num_time_slices
            }

    # In the context of HemeLB this is meaningless but a file
    # needs to be created to keep Luigi happy. As the constant pressure
    # initial condition does not need an input file any old file will do.
    def write_icond_input(self, input_target):
        self.write_input(input_target, self.initial_conditions)
        
    def write_input(self, input_target, ic):
        sub_ival = self.sub_times[input_target.i]
        res = {
            Resolution.Fine: 'fine',
            Resolution.Coarse: 'coarse'
            }[input_target.res]
            
        state = {
            'params': self.sim_params,
            'initial': ic,
            'time': sub_ival.to_dict(),
            'solver': input_target.res.name.lower()
            }
        with input_target.open('w') as f:
            yaml.dump(state, stream=f)

    def final_task(self):
        return y.producer(Resolution.Fine,
                              self.num_time_slices,
                              self.num_parareal_iters)
    
    pass


########################################################################
# Targets
########################################################################

class ParaRealTarget(luigi.LocalTarget):
    '''A local file pattern, parameterised on resolution, i, k.

    Subclass must set _pattern before this constructor is called.
    '''
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

    def __str__(self):
        return '{}[res={res.name}, i={i:d}, k={k:d}]'.format(type(self).__name__,
                                                                 **vars(self))
    
    @classmethod
    @abstractmethod
    def producer(cls, res, i, k):
        '''Targets must provide a class method producer(cls, res, i, k)
        which will return the luigi.Task (subclass) instance that will
        generate cls(res, i, k) as their output.
        '''
        pass
    
    pass

def check_producer(prod_func):
    '''Sanity check on the producer class methods that the task does
    have an output of the required type.

    Zero runtime cost in non-debug mode.
    '''
    
    if __debug__:
        @wraps(prod_func)
        def wrapper(cls, res, i, k):
            assert i >= 0, 'Time slice index must be non-negative'
            assert k >= 0, 'Parareal iteration index must be non-negative'
            assert k <= i, 'Parareal iteration must not be greater than time slice'
            
            tsk = prod_func(cls, res, i, k)
            assert tsk.output() == cls(res, i, k)
                
            return tsk
        return wrapper
    else:
        return prod_func

def logrun(run_func):
    '''Helper to log run events consistently'''
    logger = logging.getLogger('parareal')
    @wraps(run_func)
    def wrapper(self):
        target = self.output()
        logger.info('%s -> %s', type(self).__name__, target)
        return run_func(self)
    return wrapper

class Input(ParaRealTarget):
    '''Represent an input file, ready to be run.
    '''
    _pattern ='{resolution}_{i:04d}_{k:02d}/checkpoint.xtr'

    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        if i == 0:
            return InitialConditionMaker(res)
        else:
            return InputMaker(res=res, i=i, k=k)
    pass

class Trajectory(ParaRealTarget):
    '''Represent the results of running a solver.
    '''
    _pattern = '{resolution}_{i:04d}_{k:02d}/results/Extracted/checkpoint.xtr'
    
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        return Solver(res=res, i=i, k=k)

    def get_last(self):
        pass
    
    pass

class State(ParaRealTarget):
    '''A state is the state vector of the system.
    
    For now this is just a yaml file containing t, x, v (floats)
    '''

    def load(self):
        with self.open() as reader:
            return yaml.load(reader)

    def save(self, obj):
        with self.open('w') as writer:
            yaml.dump(obj, stream=writer)

    pass

class y(State):
    '''Represent a state vector that is the value of y at a given time &
    parareal iteration.
    '''
    
    _pattern = 'state/y_{resolution}_{i:04d}_{k:02d}_checkpoint.xtr'
            
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        assert i > 0, 'Should probably never request y[i=0] explicitly'
        
        if k == 0 or k == i:
            # y[i, k=0] = G(y[i-1, k=0])
            # y[i,k=i] = F(y(i=i-1, k=i-1))
            return Assignment(res=res, i=i, k=k)
        else:
            return ParaUpdate(res=res, i=i, k=k)
    
    pass

class Gy(State):
    '''A state vector that is the result of applying the coarse time
    evolution operator to y_i^k.
    '''
    
    _pattern = 'state/Gy_{resolution}_{i:04d}_{k:02d}_checkpoint.xtr'

    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src_cls = {
            Resolution.Fine: SolutionRefinement,
            Resolution.Coarse: SolutionExtractor
            }[res]
        return src_cls(res=Resolution.Coarse, i=i, k=k)
        
    pass

class DeltaGy(State):
    '''A state vector representing:
    G(y_i^k) - G(y_i^{k-1})
    '''
    
    _pattern = 'state/DGy_{resolution}_{i:04d}_{k:02d}_checkpoint.xtr'
    
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src_cls = {
            Resolution.Fine: DiffRefinement,
            Resolution.Coarse: CoarseDiffer
            }[res]
        return src_cls( i=i, k=k)
    pass

class Fy(State):
    '''A state vector that is the result of applying the fine time
    evolution operator to y_i^k.
    '''
    
    _pattern = 'state/Fy_{resolution}_{i:04d}_{k:02d}_checkpoint.xtr'
    
    @classmethod
    @check_producer
    def producer(cls, res, i, k):
        src_cls = {
            Resolution.Fine: lambda i,k: SolutionExtractor(res=Resolution.Fine, i=i, k=k),
            Resolution.Coarse: lambda i,k: SolutionCoarsening(i=i, k=k)
            }[res]
        return src_cls(i, k)
    
    pass

######################################################################
# Tasks
######################################################################

class InitialConditionMaker(luigi.ExternalTask):
    '''Produce the initial conditions.'''
    
    # Which resolution 
    res = luigi.EnumParameter(enum=Resolution)
    
    def output(self):
        return Input(self.res, 0, 0)
    @logrun
    def run(self):
        print('InitialConditionMaker:')
        print('output:')
        print(self.output().path)
        parareal().write_icond_input(self.output())
        if self.res == Resolution.Fine:
            print('Trying to link fine_initial_config.xml to fine_0000_00/config.xml')
            os.link('fine_initial_config.xml', 'fine_0000_00/config.xml')
        else:
            print('Trying to link coarse_initial_config.xml to coarse_0000_00/config.xml')
            os.link('coarse_initial_config.xml', 'coarse_0000_00/config.xml')
        return
    pass


class Op(luigi.Task):
    '''Base class for parareal operations parameterised by time slice,
    parareal iteration and resolution.

    All are required.
    '''
    
    # time slice ID
    i = luigi.IntParameter()
    # parareal ID
    k = luigi.IntParameter()
    # Which resolution 
    res = luigi.EnumParameter(enum=Resolution)
    
    pass

class InputMaker(Op):
    '''Produce all things necessary to run a time evolution from a given
    y state.
    '''
    def output(self):
        return Input(self.res, self.i, self.k)
    
    def requires(self):
        return y.producer(self.res, self.i, self.k)
    @logrun
    def run(self):
        print(self)
        print('input:')
        print(self.input().path)
        print('output:')
        print(self.output().path)
        if (self.res == Resolution.Coarse):
            dirStr = '{}_{:04d}_{:02d}'.format('coarse', self.i, self.k)
            os.makedirs(dirStr)
            print('Trying to link coarse_restart_config.xml to ' + dirStr + '/config.xml')
            os.link('coarse_restart_config.xml', dirStr + '/config.xml')
        else:
            dirStr = '{}_{:04d}_{:02d}'.format('fine', self.i, self.k)
            os.makedirs(dirStr)
            print('Trying to link fine_restart_config.xml to ' + dirStr + '/config.xml')
            os.link('fine_restart_config.xml', dirStr + '/config.xml')
        print('dirStr = ' + dirStr)
        inStr = self.input().path
        inStr1 = self.input().path[:-3] + 'off'
        outStr = self.output().path
        outStr1 = self.output().path[:-3] + 'off'
        print('InputMaker:Trying to link ' + inStr + ' to ' + outStr)
        os.link(inStr, outStr)
        os.link(inStr1, outStr1)
        return
    pass

class Solver(Op):
    '''Do one iteration of the basic solvers.
    
    Advances time from i -> i+1
    '''
    
    def output(self):
        return Trajectory(self.res, self.i, self.k)
    
    def requires(self):
        return Input.producer(self.res, self.i, self.k)
    @logrun
    def run(self):
        sys.stdout.write('Run:Solver: ')
        print(self)
        sys.stdout.write('Run:Input path: ' + self.input().path + '\n')
        sys.stdout.write('Run:Output path: ' + self.output().path + '\n')
        if (self.res == Resolution.Coarse):
            dirStr = '{}_{:04d}_{:02d}'.format('coarse', self.i, self.k)
            numProcs = numCoarseProcs
        else:
            dirStr = '{}_{:04d}_{:02d}'.format('fine', self.i, self.k)
            numProcs = numFineProcs
        with cd(dirStr):
            run_hemelb(numProcs)
            print('HemeLB output: ' + dirStr + '/results/Extracted/checkpoint.xtr')
        return
    
    pass

class SolutionExtractor(Op):
    '''Given a trajectory, pull out the state to a separate target, G(y)
    or F(y) depending on our resolution.
    '''
    def output(self):
        outcls = {Resolution.Coarse: Gy, Resolution.Fine: Fy}[self.res]
        return outcls(self.res, self.i, self.k)
    
    def requires(self):
        return Trajectory.producer(self.res, self.i, self.k)
    @logrun
    def run(self):
        # In general the output from a simulation will have data for a number of time steps.
        # The solution extractor extracts the state at the last time step.
        # In the current HemeLB implementation there must be only one state in the trajectory
        # so no extraction is required but in general 'get_last' in the class Trajectory
        # should have a non-trivial implementation and should be called here.
        print('SolutionExtractor:')
        print(self.input().path)
        print(self.output().path)
        inStr = self.input().path
        inStr1 = self.input().path[:-3] + 'off'
        outStr = self.output().path
        outStr1 = self.output().path[:-3] + 'off'
        print('SolnEx:Trying to link ' + inStr + ' to ' + outStr)
        os.link(inStr, outStr)
        os.link(inStr1, outStr1)
        return
    pass

class Assignment(Op):
    '''Since we have assignments like y[i, k=0] = G(y[i-1, k=0]) in the
    description of the parareal algorithm, we can either copy the state
    or use filesystem links. Choosing to use hard links to efficiency
    here.
    '''
    def requires(self):
        if self.k == 0:
            # y[i, k=0] = G(y[i-1, k=0])
            if self.res == Resolution.Coarse:
                return SolutionExtractor(res=Resolution.Coarse, i=self.i-1, k=0)
            else:
                return SolutionRefinement(i=self.i-1, k=0)
        
        elif self.k == self.i:
            # y[i,k=i] = F(y(i=i-1, k=i-1))
            if self.res == Resolution.Fine:
                return SolutionExtractor(res=Resolution.Fine, i=self.i-1, k=self.k-1)
            else:
                return SolutionCoarsening(i=self.i-1, k=self.k-1)
        else:
            raise ValueError("should never get here")
        
    def output(self):
        return y(self.res, self.i, self.k)
    @logrun
    def run(self):
        print('Assignment:')
        print('input:')
        print(self.input().path)
        print('output:')
        print(self.output().path)
        inStr = self.input().path
        inStr1 = self.input().path[:-3] + 'off'
        outStr = self.output().path
        outStr1 = self.output().path[:-3] + 'off'
        print('SolnEx:Trying to link ' + inStr + ' to ' + outStr)
        os.link(inStr, outStr)
        os.link(inStr1, outStr1)
        return

class Coarsening(luigi.Task):
    '''Base class for coarsening operations.

    Subclasses must define state_cls attribute
    
    No need for a resolution parameter as the input will always be Fine
    and the output Coarse.
    '''
    i = luigi.IntParameter()
    k = luigi.IntParameter()
    
    def requires(self):
        return self.state_cls.producer(Resolution.Fine, self.i, self.k)
    
    def output(self):
        return self.state_cls(Resolution.Coarse, self.i, self.k)
    @logrun
    def run(self):
        print('Coarsening:')
        print(self)
        print('input:')
        print(self.input().path)
        print('output:')
        print(self.output().path)
        fineCheckpointFileName = self.input().path
        coarseCheckpointFileName = self.output().path
        if os.path.isfile(fineCheckpointFileName):
            coarsen.createCoarseCheckpointFile(fineCheckpointFileName,
                                               coarseRankFileName,
                                               coarseCheckpointFileName)
        else:
            sys.stderr.write('ERROR: '+ fineCheckpointFileName + ' not found.\n')
        return
    pass

class SolutionCoarsening(Coarsening):
    state_cls = Fy
    pass

class Refinement(luigi.Task):
    '''Base class for refinement operations.

    Subclasses must define state_cls attribute
    
    No need for a resolution parameter as the input will always be Coarse
    and the output Fine.
    '''
    i = luigi.IntParameter()
    k = luigi.IntParameter()
    
    def requires(self):
        return self.state_cls.producer(Resolution.Coarse, self.i, self.k)
    
    def output(self):
        return self.state_cls(Resolution.Fine, self.i, self.k)
    @logrun
    def run(self):
        print('Refinement:')
        print(self)
        print('input:')
        print(self.input().path)
        print('output:')
        print(self.output().path)
        coarseCheckpointFileName = self.input().path
        fineCheckpointFileName = self.output().path
        if os.path.isfile(coarseCheckpointFileName):
            refine.createFineCheckpointFile(coarseCheckpointFileName,
                                            fineRankFileName,
                                            fineCheckpointFileName)
        else:
            sys.stderr.write('ERROR: '+ coarseCheckpointFileName + ' not found.\n')
        return
    
class SolutionRefinement(Refinement):
    '''Refine a solution.'''
    state_cls = Gy
    pass

class DiffRefinement(Refinement):
    '''Refine the difference between two coarse solutions.'''
    state_cls = DeltaGy
    pass

class CoarseDiffer(luigi.Task):
    '''Perform the G(y_i^k) - G(y_i^{k-1}) operation needed in the
    parareal update.

    Resolution will always be Coarse.
    '''
    i = luigi.IntParameter()
    k = luigi.IntParameter()
    
    def requires(self):
        return [
            Gy.producer(Resolution.Coarse, self.i, self.k),
            Gy.producer(Resolution.Coarse, self.i, self.k-1)
            ]
    
    def output(self):
        return DeltaGy(Resolution.Coarse, self.i, self.k)
    @logrun
    def run(self):
        print('CoarseDiffer:')
        print('inputs:')
        for i in self.input():
            print(i.path)
        print('output:')
        print(self.output().path)
        difference.createDifferenceDistributionFile(self.input()[0].path,
                                                self.input()[1].path,
                                                self.output().path)
        return
    pass

class ParaUpdate(Op):
    '''Do the parareal update:

    y_i^k = G(y_{i-1}^k) + (F(y_{i-1}^{k-1}) - G(y_{i-1}^{k-1}))

    But noting that the subtraction of the G(y)'s can be more efficiently
    done at the coarse resolution before refinement.
    '''
    def output(self):
        return y(self.res, self.i, self.k)
    
    def requires(self):
        return [
            DeltaGy.producer(self.res, self.i-1, self.k),
            Fy.producer(self.res, self.i-1, self.k-1)
            ]
    @logrun
    def run(self):
        print('ParaUpdate:')
        print('inputs:')
        for i in self.input():
            print(i.path)
        print('output:')
        print(self.output().path)
        sum.createSummedDistributionFile(self.input()[0].path,
                                            self.input()[1].path,
                                            self.output().path)
        return
    pass



if __name__ == '__main__':
    num_workers = int(sys.argv[1])
    os.makedirs('state')
    tsk = y.producer(Resolution.Fine, 4, 2)
    # #tsk = Solver(res=Resolution.Coarse, i=3, k=0)
    luigi.build([tsk], workers=num_workers, local_scheduler=True)
