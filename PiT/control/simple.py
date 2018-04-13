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
            return yaml.dump(obj, stream=writer)
        
    # @property
    # def val(self):
    #     return self._db.get((self.i, self.k), None)

    # @val.setter
    # def val(self, v):
    #     self._db[(self.i, self.k)] = v
    #     return

    pass

class y(State):
    _db = {}
    prefix = 'y'
    @staticmethod
    def producer(i, k):
        if i == 0:
            return InitialCondition()
        else:
            if k == 0:
                return Coarse(i=i, k=k)
            else:
                return ParaUpdate(i=i, k=k)
    def __str__(self):
        return "y[i={}, k={}]".format(self.i, self.k)
    pass

class Gy(State):
    _db = {}
    prefix = 'Gy'
    @staticmethod
    def producer(i, k):
        return Coarse(i=i, k=k)
    def __str__(self):
        return "G(y[i={}, k={}])".format(self.i, self.k)
    pass

class Fy(State):
    _db = {}
    prefix = 'Fy'
    @staticmethod
    def producer(i, k):
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
        self.output().save(yik.load())
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
        self.output().save(yik.load())
    pass

class ParaUpdate(Op):
    def output(self):
        return y(self.i, self.k)

    def requires(self):
        return [Gy.producer(self.i-1, self.k), Fy.producer(self.i-1, self.k-1), Gy.producer(self.i-1, self.k-1)]
        
    def run(self):
        g_cur, f_last, g_last = self.input()
        self.output().save(g_cur.load() + f_last.load() - g_last.load())
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
        
    num_time_slices = 12
    num_parareal_iters = 3

    ic = y(0,0)
    ic.save(1.0)
        
    final_task = y.producer(num_time_slices, num_parareal_iters)
    luigi.build([final_task])

    logfile.close()
