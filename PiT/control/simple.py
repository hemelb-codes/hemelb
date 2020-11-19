# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import luigi
import yaml

class State(luigi.LocalTarget):
    def __init__(self, i, k):
        self.i = i
        self.k = k
        path = '{}_{:04d}_{:02d}'.format(self.prefix, i, k)
        super(State, self).__init__(path)
        return
    
    def load(self):
        with self.open() as reader:
            return yaml.load(reader)

    def save(self, obj):
        with self.open('w') as writer:
            yaml.dump(obj, stream=writer)
        
    pass

class y(State):
    prefix = 'y'
    @staticmethod
    def producer(i, k):
        assert i >= 0
        assert k >= 0
        assert k <= i
        
        if i == 0:
            return InitialCondition()
        else:
            if k == 0:
                return Coarse(i=i-1, k=k)
            elif k == i:
                return Fine(i=i-1, k=k-1)
            else:
                return ParaUpdate(i=i, k=k)
    def __str__(self):
        return "y[i={}, k={}]".format(self.i, self.k)
    pass

class Gy(State):
    prefix = 'Gy'
    @staticmethod
    def producer(i, k):
        if i == 0:
            k = 0
        return Coarse(i=i, k=k)
    def __str__(self):
        return "G(y[i={}, k={}])".format(self.i, self.k)
    pass

class Fy(State):
    prefix = 'Fy'
    @staticmethod
    def producer(i, k):
        if i == 0:
            k = 0
        return Fine(i=i, k=k)
    def __str__(self):
        return "F(y[i={}, k={}])".format(self.i, self.k)

    pass

class Op(luigi.Task):
    # These identify the state after execution
    
    # time slice ID
    i = luigi.IntParameter()
    # parareal ID
    k = luigi.IntParameter()

class Coarse(Op):
    def output(self):
        return Gy(self.i, self.k)
    
    def requires(self):
        return y.producer(self.i, self.k)
    
    def run(self):
        log(self.output())
        yik = self.input()
        data = yik.load()
        assert data['t'] == self.i
        data['t'] += 1
        self.output().save(data)
        return
    pass

class Fine(Op):
    def output(self):
        return Fy(self.i, self.k)
    
    def requires(self):
        return y.producer(self.i, self.k)
    
    def run(self):
        log(self.output())
        yik = self.input()
        data = yik.load()
        assert data['t'] == self.i
        data['t'] += 1
        self.output().save(data)
        return
    pass

class ParaUpdate(Op):
    def output(self):
        return y(self.i, self.k)

    def requires(self):
        return [Gy.producer(self.i-1, self.k), Fy.producer(self.i-1, self.k-1), Gy.producer(self.i-1, self.k-1)]
        
    def run(self):
        inputs = self.input()
        idata = [i.load() for i in inputs]        
        uniq_t = set(i['t'] for i in idata)
        assert len(uniq_t) == 1
        assert uniq_t.pop() == self.i

        g_cur, f_last, g_last = [i['y'] for i in idata]
        ans = g_cur + f_last - g_last
        
        self.output().save({'t': self.i, 'y': ans})
        
        g_cur, f_last, g_last = inputs
        log("{} = {} + {} - {}".format(self.output(), g_cur, f_last, g_last))
    pass

class InitialCondition(luigi.ExternalTask):
    def output(self):
        return y(0,0)
    pass

if __name__ == "__main__":
    logfile = open('run.log', 'w')
    def log(*things):
        logfile.write(' '.join(str(t) for t in things) + '\n')
        
    num_time_slices = 4
    num_parareal_iters = 2

    ic = y(0,0)
    ic.save({'t': 0, 'y': 1.0})
        
    final_task = y.producer(num_time_slices, num_parareal_iters)
    luigi.build([final_task], local_scheduler=True)
    logfile.close()
