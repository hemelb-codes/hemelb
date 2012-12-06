import pytest

import itertools
import math
import time

from HemeLbSetupTool.Util.Observer import Observable, NotifyOptions, ObservableList

class Circle(Observable):
    def __init__(self, r):
        self.radius = r
        self.AddDependency('diameter', 'radius')
        self.AddDependency('perimeter', 'diameter')
        self.AddDependency('area', 'radius')
        return
    
    @property
    def diameter(self):
        return 2*self.radius
    @property
    def perimeter(self):
        return math.pi * self.diameter
    @property
    def area(self):
        return math.pi * self.radius**2
    pass

class Something(Observable):
    def __init__(self, r):
        self.circle = Circle(r)
        return
    pass

class Observer(object):
    """Helper - will record when it was called as precisely as the system allows"""
    def __init__(self):
        self.CallTime = None
        self.Change = None
        return
    def __call__(self, change):
        self.CallTime = time.clock()
        self.Change = change
        time.sleep(0.01)
        return
    def Reset(self):
        self.CallTime = None
        self.Change = None
    pass

class TestObservable:
    def test_circle(self):
        c = Circle(5.)
        assert c.diameter == 10.
        assert c.perimeter == 10 * math.pi
        assert c.area == 25 * math.pi
    
    def test_observation(self):
        # Check it works at all
        ob = Observer()
        
        c = Circle(5.)
        c.AddObserver('radius', ob)
        c.radius = 10.
        
        assert ob.CallTime is not None
        
        # Check the ordering of PRE and POST change is correct
        preOb = Observer()
        postOb = Observer()
        c.AddObserver('radius', preOb, options=NotifyOptions(BEFORE_CHANGE=True,
                                                             AFTER_CHANGE=False))
        c.AddObserver('radius', postOb, options=NotifyOptions(BEFORE_CHANGE=False,
                                                              AFTER_CHANGE=True))
        
        c.radius = 5.
        
        assert preOb.CallTime < postOb.CallTime
        return
    
    def test_dependency(self):
        c = Circle(5.)
        ob = Observer()
        c.AddObserver('perimeter', ob)
        c.radius = 10.
        assert ob.CallTime is not None
        return
    
    def test_keypath(self):
        # Test the basic working of complex observation keypaths 
        obj = Something(7.)
        ob = Observer()
        obj.AddObserver('circle.diameter', ob)
        obj.circle.radius = 5.
        
        assert ob.CallTime is not None
        
        # Test that this correctly unsubscribes the old object and subscribes the new one
        ob.Reset()
        oldCirc = obj.circle
        newCirc = Circle(1)
        
        obj.circle = newCirc
        assert ob.CallTime is not None
        
        ob.Reset()
        oldCirc.radius = 77
        assert ob.CallTime is None
        
        ob.Reset()
        obj.circle.radius = 123
        assert ob.CallTime is not None
        
        return
    
    pass

class TestObservableList:
    def test_listyness(self):
        old = range(10)
        new = ObservableList(old)
        
        assert old == new
        
        # Make them differ in first elem
        old[0] = 17
        assert old != new
        
        assert old[3] == new[3]
        
        assert old[3:9:2] == new[3:9:2]
        
        with pytest.raises(IndexError):
            new[17]
        
        # Delete 1st elem, should be the same again
        del new[0]
        del old[0]
        assert old == new
        
        # Check len
        assert len(new) == 9
        assert new.length == 9
        
        new.insert(0, 'foo')
        old.insert(0, 'foo')
        assert old == new
        
        return
    
    @staticmethod
    def MakeListObservers(lst):
        pre = NotifyOptions(BEFORE_CHANGE=True,
                            AFTER_CHANGE=False)
        post = NotifyOptions(BEFORE_CHANGE=False,
                             AFTER_CHANGE=True)
        
        preInsOb = Observer()
        lst.AddObserver('@INSERTION', preInsOb, pre)
        postInsOb = Observer()
        lst.AddObserver('@INSERTION', postInsOb, post)
        
        preRemOb = Observer()
        lst.AddObserver('@REMOVAL', preRemOb, pre)
        postRemOb = Observer()
        lst.AddObserver('@REMOVAL', postRemOb, post)
        
        preRepOb = Observer()
        lst.AddObserver('@REPLACEMENT', preRepOb, pre)
        postRepOb = Observer()
        lst.AddObserver('@REPLACEMENT', postRepOb, post)
        
        preLenOb = Observer()
        lst.AddObserver('length', preLenOb, pre)
        postLenOb = Observer()
        lst.AddObserver('length', postLenOb, post)
        return ((preInsOb, postInsOb),
                (preRemOb, postRemOb),
                (preRepOb, postRepOb),
                (preLenOb, postLenOb))
    
    
                    
    def test_observation(self):
        def ResetAll(obs):
            for obPair in obs:
                for ob in obPair:
                    ob.Reset()
            return
    
        def CheckObservers(observers, pattern):
            for obPair, shouldChange in itertools.izip(observers, pattern):
                if shouldChange:
                    for ob in obPair:
                        assert ob.CallTime is not None
                    assert obPair[0].CallTime < obPair[1].CallTime
                else:
                    for ob in obPair:
                        assert ob.CallTime is None
            return
        
        lst = ObservableList()
        
        observers = self.MakeListObservers(lst)
        
        # Add an elem
        lst.append(Circle(4))
        CheckObservers(observers, (True, False, False, True))
        
        # Replace one
        ResetAll(observers)
        lst[0] = 18
        CheckObservers(observers, (False, False, True, False))
        
        # Delete one
        ResetAll(observers)
        del lst[0]
        CheckObservers(observers, (False, True, False, True))
        
        return
    
    def test_copy(self):
        import copy
        old = ObservableList(range(3))
        new = copy.copy(old)
        
        new.append('wibble')
        
        with pytest.raises(IndexError):
            old[3]
        
        
        
    pass
#olist = ObservableList()
#print 'Add b as observer of olist.@INSERTION'
#olist.AddObserver('@INSERTION', b.change)
#print 'Append to olist'
#olist.append(6)
