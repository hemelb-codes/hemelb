# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

"""
To demonstrate how HemeLB may be called from Python in the context
of the ParaReal algorithm.
"""

import luigi
import sys
import yaml
import time
from contextlib import contextmanager
import os
import subprocess
from shutil import copyfile
from string import Template

import coarsen
import refine
import correction

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
        sys.stderr.write(' ' + item)
    sys.stderr.write('\n')
    p = subprocess.Popen(runcommand, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdoutdata, stderrdata = p.communicate()
    sys.stdout.write(stdoutdata)
    sys.stderr.write(stderrdata)
    pass

useHemeLB = True

numFineProcs = 4
numCoarseProcs = 3

coarseRankFileName = 'initialisation/coarseMpirank.xtr'
fineRankFileName = 'initialisation/fineMpirank.xtr'

#*******************************************************************************
# Targets
#*******************************************************************************

class State(luigi.LocalTarget):
    def __init__(self, i, k):
        self.i = i
        self.k = k
        self.dir = '{}_{:04d}_{:02d}'.format(self.prefix, i, k)
        path = self.dir + '/file'
        super(State, self).__init__(path)
        return
    
    def load(self):
        with self.open() as reader:
            return yaml.load(reader)

    def save(self, obj):
        with self.open('w') as writer:
            yaml.dump(obj, stream=writer)
    pass

class yC(State):
    prefix = 'yC'
    @staticmethod
    def producer(i, k):
        assert k >= 0

        # Here C is the coarsening operator.
        if k == 0:
            if i == 0:
                # yC(0, 0) = C(y(0, 0))
                return C(0, 0)
            else:
                # yC(i, 0) = Gy(i-1, 0) with input yC(i-1, 0)
                return Gy(i-1, 0)
        else:
            if i == 0:
                # yC(0, k) = C(y(0, 0))
                return C(0, 0)
            else:
                # yC(i, k) =  Gy(i-1, k) with input C(y(i-1, k))
                return Gy(i-1, k)

    def __str__(self):
        return self.dir
    pass

class yF(State):
    prefix = 'yF'
    @staticmethod
    def producer(i, k):

        if i == 0:
            return InitialCondition()
        else:
            # yF(i, k) =  Fy(i-1, k) with input y(i-1, k)
            return Fy(i-1, k)

    def __str__(self):
        return self.dir
    pass

class y(State):
    prefix = 'y'
    @staticmethod
    def producer(i, k):

        if k == 0:
            if i == 0:
                # y(0, 0) = initial condition
                return InitialCondition()
            else:
                # y(i, 0) = R(yC(i, 0))
                return R(i, 0)
        elif i == 0:
            # y(0, k) = y(0, 0)
            return InitialCondition()
        elif i < k:
            return yF(i, k-1)
        else:
            # y(i, k) = RyC(i, k) + yF(i, k-1) - RyC(i, k-1)
            return StateCorrection(i, k)

    def __str__(self):
        return self.dir
    pass

class RyC(State):
    prefix = 'RyC'
    @staticmethod
    def producer(i, k):
        # Here R is the refinement operator.
        # RyC(i, k) = R(yC(i, k))
        return R(i, k)

    def __str__(self):
        return self.dir
    pass

#*******************************************************************************
# Tasks
# (Non-external) tasks must implement:
#   run
#   output (returns targets only)
#   requires (specifies tasks only)
#*******************************************************************************

class InitialCondition(luigi.ExternalTask):
    def output(self):
        return y(0, 0)
    pass

class Op(luigi.Task):
    # Time slice ID.
    i = luigi.IntParameter()
    # ParaReal iteration index.
    k = luigi.IntParameter()

class Gy(Op):
    # Coarse integrator.
    def output(self):
        #assert self.i > 0
        return yC(self.i+1, self.k)

    def requires(self):
        return yC.producer(self.i, self.k)

    def run(self):
        # yC(i, 0) = Gy(i-1, 0) with input yC(i-1, 0)
        # yC(i, k) =  Gy(i-1, k) with input C(y(i-1, k))
        CyCik = self.input()
        data = CyCik.load()
        cwd = os.getcwd()
        sys.stderr.write("Gy, cwd: " + cwd + "\n")
        inputStr = CyCik.__str__() + '/results/Extracted/checkpoint.xtr'
        sys.stderr.write("Gy, input: " + inputStr + "\n")
        if not os.path.isfile(inputStr):
            sys.stderr.write('Gy, ERROR: ' + inputStr + ' not found.\n')
        self.output().save(data)
        configTemplateFile = open('coarse_restart_template.xml', 'r')
        configTemplateStr = configTemplateFile.read()
        configTemplateFile.close()
        configTemplate = Template(configTemplateStr)
        inputStr = os.path.abspath(inputStr)
        configStr = configTemplate.substitute(checkpoint_in=inputStr)
        dirStr = self.output().__str__()
        sys.stderr.write('Gy, spam: ' + dirStr + '\n')
        if useHemeLB:
            with cd(dirStr):
                configFile = open('config.xml', 'w')
                configFile.write(configStr)
                configFile.close()
                run_hemelb(numCoarseProcs)
        else:  # Fake the production of the output.
            longStr = dirStr + '/results/Extracted/checkpoint'
            if not os.path.exists(longStr + '.xtr'):
                cwd = os.getcwd()
                sys.stderr.write('Gy, cwd: ' + cwd + ', creating ' + longStr + '\n')
                os.makedirs(dirStr + '/results/Extracted')
                copyfile('checkpoint_coarse.xtr', longStr + '.xtr')
                copyfile('checkpoint_coarse.off', longStr + '.off')
        return
    pass

class Fy(Op):
    # Fine integrator.
    def output(self):
        #assert self.i > 0
        return yF(self.i+1, self.k)

    def requires(self):
        return y.producer(self.i, self.k)

    def run(self):
        # yF(i, k) =  Fy(i-1, k) with input y(i-1, k)
        yik = self.input()
        data = yik.load()
        cwd = os.getcwd()
        sys.stderr.write("Gy, cwd: " + cwd + "\n")
        inputStr = yik.__str__() + '/results/Extracted/checkpoint.xtr'
        sys.stderr.write("Fy, input: " + inputStr + "\n")
        if not os.path.isfile(inputStr):
            sys.stderr.write('Fy, ERROR: '+ inputStr + ' not found.\n')
        self.output().save(data)
        configTemplateFile = open('fine_restart_template.xml', 'r')
        configTemplateStr = configTemplateFile.read()
        configTemplateFile.close()
        configTemplate = Template(configTemplateStr)
        inputStr = os.path.abspath(inputStr)
        configStr = configTemplate.substitute(checkpoint_in=inputStr)
        dirStr = self.output().__str__()
        sys.stderr.write('Fy, spam: ' + dirStr + '\n')
        if useHemeLB:
            with cd(dirStr):
                configFile = open('config.xml','w')
                configFile.write(configStr)
                configFile.close()
                run_hemelb(numFineProcs)
        else: # Fake the production of the output.
            longStr = dirStr + '/results/Extracted/checkpoint'
            if not os.path.exists(longStr + '.xtr'):
                cwd = os.getcwd()
                sys.stderr.write('Fy, cwd: ' + cwd + ', creating ' + longStr + '\n')
                os.makedirs(dirStr + '/results/Extracted')
                sys.stderr.write('Fy: creating ' + longStr + '\n')
                copyfile('checkpoint_fine.xtr', longStr + '.xtr')
                copyfile('checkpoint_fine.off', longStr + '.off')
        return
    pass

class StateCorrection(Op):
    def output(self):
        assert self.i > 0
        return y(self.i, self.k)

    def requires(self):
        return [RyC.producer(self.i, self.k), yF.producer(self.i, self.k-1), RyC.producer(self.i, self.k-1)]

    def run(self):
        # y(i, k) = RyC(i, k) + yF(i, k-1) - RyC(i, k-1)
        inputs = self.input()
        checkpointFileNames = []
        for i in inputs:
            fileName = i.__str__() + '/results/Extracted/checkpoint.xtr'
            sys.stderr.write("SC, input: " + fileName + "\n")
            checkpointFileNames.append(fileName)
            if not os.path.isfile(fileName):
                sys.stderr.write('SC, ERROR: '+ fileName + ' not found.\n')

        data = {'state': 'y', 't': self.i, 'y': 1}
        self.output().save(data)
        
        correctedCheckpointFileName = self.output().__str__() + '/results/Extracted/checkpoint.xtr'
        correction.createCorrectedDistributionFile(checkpointFileNames[0],
                                                   checkpointFileNames[1],
                                                   checkpointFileNames[2],
                                                   correctedCheckpointFileName)
pass

class R(Op):
    # R is the refinement operator.
    def output(self):
        return y(self.i, self.k)

    def requires(self):
        return yC.producer(self.i, self.k)

    def run(self):
        # RyC(i, k) = R(yC(i, k))
        yCik = self.input()
        data = yCik.load()
        self.output().save(data)
        coarseCheckpointFileName = yCik.__str__() + '/results/Extracted/checkpoint.xtr'
        sys.stderr.write("R, input: " + coarseCheckpointFileName + "\n")
        fineCheckpointFileName = self.output().__str__() + '/results/Extracted/checkpoint.xtr'
        if os.path.isfile(coarseCheckpointFileName):
            refine.createFineCheckpointFile(coarseCheckpointFileName,
                                            fineRankFileName,
                                            fineCheckpointFileName)
        else:
            sys.stderr.write('R, ERROR: '+ coarseCheckpointFileName + ' not found.\n')
pass

class C(Op):
    # C is the coarsening operator.
    def output(self):
        return yC(self.i, self.k)

    def requires(self):
        return y.producer(self.i, self.k)

    def run(self):
        # yC(i, k) = C(y(i, k))
        yik = self.input()
        data = yik.load()
        self.output().save(data)
        fineCheckpointFileName = yik.__str__() + '/results/Extracted/checkpoint.xtr'
        sys.stderr.write("C, input: " + fineCheckpointFileName + "\n")
        coarseCheckpointFileName = self.output().__str__() + '/results/Extracted/checkpoint.xtr'
        if os.path.isfile(fineCheckpointFileName):
            coarsen.createCoarseCheckpointFile(fineCheckpointFileName,
                                               coarseRankFileName,
                                               coarseCheckpointFileName)
        else:
            sys.stderr.write('C, ERROR: '+ fineCheckpointFileName + ' not found.\n')
pass

#*******************************************************************************
# Entry Point
#*******************************************************************************

if __name__ == "__main__":
    logfile = open('run.log', 'w')
    def log(*things):
        logfile.write(' '.join(str(t) for t in things) + '\n')

    print "num_time_slices =", sys.argv[1]
    print "num_parareal_iters =", sys.argv[2]
    print "num_workers =", sys.argv[3]
        
    num_time_slices = int(sys.argv[1])
    num_parareal_iters = int(sys.argv[2])
    num_workers = int(sys.argv[3])

    ic = y(0, 0)
    ic.save({'state': 'y', 't': 0, 'y': 1.0})

    os.makedirs('y_0000_00/results/Extracted')
    copyfile('initialisation/checkpoint_fine.xtr', 'y_0000_00/results/Extracted/checkpoint.xtr')
    copyfile('initialisation/checkpoint_fine.off', 'y_0000_00/results/Extracted/checkpoint.off')
    if num_parareal_iters == 0:
        final_task = yC.producer(num_time_slices, num_parareal_iters)
    else:
        final_task = y.producer(num_time_slices, num_parareal_iters)
    luigi.build([final_task], workers=num_workers, local_scheduler=True)
    logfile.close()
