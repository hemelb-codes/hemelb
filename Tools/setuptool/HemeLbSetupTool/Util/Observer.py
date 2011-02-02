"""Implementation of the Observer pattern.

Create a subclass of Observable to use.
"""
import collections
from Enum import Enum

class NoSuch(object):
    """Special None-like object for when an attribute did not exist.
    """
    pass
NoSuch = NoSuch()


# ChangeKinds = Enum("SETTING", "INSERTION", "REMOVAL", "REPLACEMENT")
ChangeTimes = Enum("BEFORE", "AFTER")

class NotifyOptions(object):
    def __init__(self, **kwargs):
        defaults = {'BEFORE_CHANGE': False,
                    'AFTER_CHANGE': True# ,
                    # 'FOR_SETTING': True,
                    # 'FOR_INSERTION': False,
                    # 'FOR_REMOVAL': False,
                    # 'FOR_REPLACEMENT': False,
                    }

        for key, default in defaults.iteritems():
            setattr(self, key, kwargs.pop(key, default))
            continue

        if len(kwargs):
            raise TypeError("Invalid keyword argument '%s'" % key)

        return
    pass

class Change(object):
    """Describe a change.
    """

    def __init__(self, time=None, obj=None, attr=None, new=None, old=None, index=None):
        # if kind not in ChangeKinds:
        #     raise ValueError("Keyword argument 'kind' must be specified and one of %s" %
        #                      str(ChangeKinds))
        if time not in ChangeTimes:
            raise ValueError("Keyword argument 'time' must be specified and one of %s" %
                             str(ChangeTimes))
        # self.kind = kind
        self.time = time
        self.obj = obj
        self.new = new
        self.old = old
        self.index = index
        return

    pass


class Observable(object):
    """A base class/mixin to make attributes/properties of the user
    class observable. You should ensure that this base class' __new__
    is called before any other initialization. Just implementing the
    usual __init__ will take care of this.

    """
    
    def __new__(cls, *args, **kwargs):
        """Use __new__ to preempt user __init__ method.
        """
        # new = super(Observable, cls).__new__(cls, *args, **kwargs)
        new = object.__new__(cls)
        # Explicitly use object setattr to avoid our overrided one
        # that needs these attributes to have already been set!  Note
        # that we have to do the name mangling ourselves to mark them
        # as private
        object.__setattr__(new, '_Observable__preObservers', {})
        object.__setattr__(new, '_Observable__postObservers', {})
        object.__setattr__(new, '_Observable__dependencies', {})
        return new

    def __setattr__(self, name, newValue):
        """Set the attribute in the usual fashion, but notify any
        observers of this attribute.
        """
        # # First get the old value, covering the case that it doesn't yet exist.
        # try:
        #     oldValue = getattr(self, name)
        #     if newValue is oldValue:
        #         # If the objects are the same, this isn't a change so
        #         # return.
        #         return
        # except AttributeError:
        #     # This attribute hasn't previously been set, so use our
        #     # special NoSuch object
        #     oldValue = NoSuch
        #     pass
        
        self.WillChangeAttr(name)
        
        object.__setattr__(self, name, newValue)
        
        self.DidChangeAttr(name)
        return

    def __notifyList(self, oMap, attr, **changeOpts):
        try:
            # Get the set of observers
            oList = oMap[attr]
        except KeyError:
            # No observers have been set, so return
            return
        opts = dict(# kind=ChangeKinds.SETTING,
                    obj=self,
                    attr=attr)
        opts.update(changeOpts)
        # Construct a change object
        change = Change(**opts)
        
        for cb in oList:
            # Call each observer callback with the change object
            cb(change)
            continue
        return

    def WillChangeAttr(self, name, **changeOpts):
        """Tell our observes that we're about to change.
        """
        changeOpts['time'] = ChangeTimes.BEFORE
        self.__notifyList(self.__preObservers, name, **changeOpts)
        return
    
    def DidChangeAttr(self, name, **changeOpts):
        """Tell our observers that we have just changed.
        """
        changeOpts['time'] = ChangeTimes.AFTER
        self.__notifyList(self.__postObservers, name, **changeOpts)
        return
    
    def AddObserver(self, attr, callback, options=NotifyOptions()):
        """Make 'callback' an observer of changes to the attribute
        'attr'. The callback must be a callable taking one argument,
        which will be an Observer.Change object.
        """
        assert isinstance(attr, str)
        for flag, obs in ((options.BEFORE_CHANGE, self.__preObservers),
                          (options.AFTER_CHANGE, self.__postObservers)):
            
            if flag:
                try:
                    attrObs = obs[attr]
                except KeyError:
                    attrObs = obs[attr] = set()
                    pass
                attrObs.add(callback)
                pass
            continue
        
        return

    def RemoveObserver(self, attr, observer):
        """Remove the observer from the list to be notified on change.
        """
        for timeObservers in (self.__preObservers, self.__postObservers):
            try:
                attrObs = timeObservers[attr]
                attrObs.remove(observer)
            except KeyError:
                pass
            continue
        
        return
    
    def AddDependency(self, ofAttr, onAttr):
        """Mark an attribute/property (ofAttr) as dependent on the
        value of another (onAttr). Thus when onAttr is changed, the
        observers of ofAttr will be notified.
        """
        try:
            depMap = self.__dependencies[ofAttr]
        except KeyError:
            depMap = self.__dependencies[ofAttr] = dict()
            pass
        
        if onAttr in depMap:
            # Trying to add a redundant dependency
            return
        depMap[onAttr] = Dependency(self, ofAttr, onAttr)
        
        self.AddObserver(
            onAttr, depMap[onAttr],
            options=NotifyOptions(BEFORE_CHANGE=True,
                                  AFTER_CHANGE=True)
            )
        return
    
    pass

class Dependency(object):
    """Helper class for notification of dependent attributes.
    """
    def __init__(self, obj, ofAttr, onAttr):
        self.obj = obj
        self.ofAttr = ofAttr
        self.onAttr = onAttr
        return

    def __call__(self, change):
        """The object is the callback itself (i.e. it's a functor).
        """
        if change.time is ChangeTimes.BEFORE:
            self.obj.WillChangeAttr(self.ofAttr)
            pass
        if change.time is ChangeTimes.AFTER:
            self.obj.DidChangeAttr(self.ofAttr)
            pass
        return
    pass

class ObservableList(Observable, collections.MutableSequence):
    """A list whose member mutations can be observed. Do this by
    observing special attributes '@INSERTION', '@REMOVAL',
    '@REPLACEMENT'.
    
    """
    def __init__(self, iterable=None):
        if iterable is None:
            self.__contents = list()
        else:
            self.__contents = list(iterable)
            pass
        return

    def __len__(self):
        return self.__contents.__len__()
    
    def __getitem__(self, index):
        return self.__contents.__getitem__(index)
    
    def __setitem__(self, index, obj):
        if index >= len(self) or index < -len(self):
            raise IndexError('ObservableList assignment index out of range')
        
        self.WillChangeAttr('@REPLACEMENT', index=index)
        ans = self.__contents.__setitem__(index, obj)
        self.DidChangeAttr('@REPLACEMENT', index=index)
    
    def __delitem__(self, index):
        if index >= len(self) or index < -len(self):
            raise IndexError('ObservableList assignment index out of range')
        
        self.WillChangeAttr('@REMOVAL', index=index)
        self.__contents.__delitem__(index)
        self.DidChangeAttr('@REMOVAL', index=index)
        return
    
    def insert(self, index, object):
        """Insert object before index.
        """
        self.WillChangeAttr('@INSERTION', index=index)
        self.__contents.insert(index, object)
        self.DidChangeAttr('@INSERTION', index=index)
        return

    def __str__(self):
        return self.__contents.__str__()
    def __repr__(self):
        return self.__contents.__repr__()
    pass

if __name__ == "__main__":
    class Observed(Observable):
        def __init__(self, x=0):
            self.x = x
            self.AddDependency('y', 'x')
            return

        @property
        def y(self):
            return self.x
        
        pass
    
    class Observer(object):
        def changex(self, change):
            print change.__dict__
            return
        pass

    print "Intantiate Observed a"
    a = Observed(7)
    b = Observer()
    print "Add b as observer of a.y"
    a.AddObserver('y', b.changex)
    print "Set a.x"
    a.x = 1
    
    olist = ObservableList()
    olist.AddObserver('@INSERTION', b.changex)
    olist.append(6)
