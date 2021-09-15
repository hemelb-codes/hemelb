# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
import pytest

import math
import time

from HemeLbSetupTool.Util.Observer import Observable, NotifyOptions, ObservableList


class Circle(Observable):
    def __init__(self, r):
        self.radius = r
        self.AddDependency("diameter", "radius")
        self.AddDependency("perimeter", "diameter")
        self.AddDependency("area", "radius")
        return

    @property
    def diameter(self):
        return 2 * self.radius

    @property
    def perimeter(self):
        return math.pi * self.diameter

    @property
    def area(self):
        return math.pi * self.radius ** 2

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
        self.CallTime = time.perf_counter()
        self.Change = change
        time.sleep(0.01)
        return

    def Reset(self):
        self.CallTime = None
        self.Change = None

    pass


def test_circle():
    c = Circle(5.0)
    assert c.diameter == 10.0
    assert c.perimeter == 10 * math.pi
    assert c.area == 25 * math.pi


def test_observation():
    # Check it works at all
    ob = Observer()

    c = Circle(5.0)
    c.AddObserver("radius", ob)
    c.radius = 10.0

    assert ob.CallTime is not None

    # Check the ordering of PRE and POST change is correct
    preOb = Observer()
    postOb = Observer()
    c.AddObserver(
        "radius", preOb, options=NotifyOptions(BEFORE_CHANGE=True, AFTER_CHANGE=False)
    )
    c.AddObserver(
        "radius", postOb, options=NotifyOptions(BEFORE_CHANGE=False, AFTER_CHANGE=True)
    )

    c.radius = 5.0

    assert preOb.CallTime < postOb.CallTime
    return


def test_dependency():
    # Create a circle and observe its perimeter
    c = Circle(5.0)
    ob = Observer()
    c.AddObserver("perimeter", ob)
    c.radius = 10.0
    assert ob.CallTime is not None

    # Remove the observer, reset
    c.RemoveObserver("perimeter", ob)
    ob.Reset()
    assert ob.CallTime is None

    # Make a change and check it was not triggered
    c.radius = 3.0
    assert ob.CallTime is None


def test_dependency_method():
    # As above, but observer is a bound method

    # Create a circle and observe its perimeter
    c = Circle(5.0)
    ob = Observer()
    c.AddObserver("perimeter", ob.__call__)
    c.radius = 10.0
    assert ob.CallTime is not None

    # Remove the observer, reset
    c.RemoveObserver("perimeter", ob.__call__)
    ob.Reset()
    assert ob.CallTime is None

    # Make a change and check it was not triggered
    c.radius = 3.0
    assert ob.CallTime is None


def test_keypath():
    # Test the basic working of complex observation keypaths

    # Create sonmething with a Circle object as a member
    obj = Something(7.0)
    ob = Observer()

    # Observe the dependend object
    obj.AddObserver("circle.diameter", ob)

    # Set by full key
    obj.SetValueForKey("circle.radius", 4.0)
    assert ob.CallTime is not None

    # Test by direct setting on circle
    ob.Reset()
    obj.circle.radius = 5.0
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


def test_listyness():
    old = list(range(10))
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

    new.insert(0, "foo")
    old.insert(0, "foo")
    assert old == new

    return


def MakeListObservers(lst):
    pre = NotifyOptions(BEFORE_CHANGE=True, AFTER_CHANGE=False)
    post = NotifyOptions(BEFORE_CHANGE=False, AFTER_CHANGE=True)

    preInsOb = Observer()
    lst.AddObserver("@INSERTION", preInsOb, pre)
    postInsOb = Observer()
    lst.AddObserver("@INSERTION", postInsOb, post)

    preRemOb = Observer()
    lst.AddObserver("@REMOVAL", preRemOb, pre)
    postRemOb = Observer()
    lst.AddObserver("@REMOVAL", postRemOb, post)

    preRepOb = Observer()
    lst.AddObserver("@REPLACEMENT", preRepOb, pre)
    postRepOb = Observer()
    lst.AddObserver("@REPLACEMENT", postRepOb, post)

    preLenOb = Observer()
    lst.AddObserver("length", preLenOb, pre)
    postLenOb = Observer()
    lst.AddObserver("length", postLenOb, post)
    return (
        (preInsOb, postInsOb),
        (preRemOb, postRemOb),
        (preRepOb, postRepOb),
        (preLenOb, postLenOb),
    )


def test_observation():
    def ResetAll(obs):
        for obPair in obs:
            for ob in obPair:
                ob.Reset()
        return

    def CheckObservers(observers, pattern):
        for obPair, shouldChange in zip(observers, pattern):
            if shouldChange:
                for ob in obPair:
                    assert ob.CallTime is not None
                assert obPair[0].CallTime < obPair[1].CallTime
            else:
                for ob in obPair:
                    assert ob.CallTime is None
        return

    lst = ObservableList()

    observers = MakeListObservers(lst)

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


def test_copy():
    import copy

    old = ObservableList(range(3))
    new = copy.copy(old)

    new.append("wibble")

    with pytest.raises(IndexError):
        old[3]
