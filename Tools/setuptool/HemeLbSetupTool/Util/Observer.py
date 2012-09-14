"""Implementation of the Observer pattern.

Create a subclass of Observable to use.
"""
import collections
from copy import copy

from .Enum import Enum

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

    def __init__(self, time=None, obj=None, key=None, new=None, old=None, index=None):
        # if kind not in ChangeKinds:
        #     raise ValueError("Keyword argument 'kind' must be specified and one of %s" %
        #                      str(ChangeKinds))
        if time not in ChangeTimes:
            raise ValueError("Keyword argument 'time' must be specified and one of %s" %
                             str(ChangeTimes))
        # self.kind = kind
        self.time = time
        self.obj = obj
        self.key = key
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
    _Args = {}
    
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

    def __setattr__(self, name, newValue, **changeOpts):
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
        
        self.WillChangeValueForKey(name, **changeOpts)
        
        object.__setattr__(self, name, newValue)
        
        self.DidChangeValueForKey(name, **changeOpts)
        return

    def __notifyList(self, oMap, key, **changeOpts):
        try:
            # Get the set of observers
            attrObservers = oMap[key]
        except KeyError:
            # No observers have been set for attr
            attrObservers = set()
            pass
        
        if '.' not in key:
            # Only for simple keys do I want to notify @ANY observers
            try:
                anyObservers = oMap['@ANY']
            except KeyError:
                # No observers set for @ANY
                anyObservers = set()
                pass
            oList = set.union(attrObservers, anyObservers)
        else:
            oList = attrObservers
            pass
        
        # No observers, quit now
        if not len(oList): return
        
        opts = dict(obj=self,
                    key=key)
        opts.update(changeOpts)
        # Construct a change object
        change = Change(**opts)
        
        for cb in oList:
            # Call each observer callback with the change object
            cb(change)
            continue
        return

    def GetValueForKey(self, keyPath):
        keyParts = keyPath.split('.', 1)
        if len(keyParts) == 1:
            return self._GetLocalValueForKey(keyPath)
        else:
            localKey, restOfKey = keyParts
            return self._GetLocalValueForKey(localKey).GetValueForKey(restOfKey)
        return
    
    def _GetLocalValueForKey(self, key):
        return getattr(self, key)
    
    def SetValueForKey(self, keyPath, value, **changeOpts):
        keyParts = keyPath.split('.', 1)
        if len(keyParts) == 1:
            self._SetLocalValueForKey(keyPath, value, **changeOpts)
        else:
            localKey, restOfKey = keyParts
            self._GetLocalValueForKey(localKey).SetValueForKey(restOfKey, value, **changeOpts)
            pass
        return
    
    def _SetLocalValueForKey(self, key, value, **changeOpts):
        self.__setattr__(key, value, **changeOpts)
        return
    
    def WillChangeValueForKey(self, key, **changeOpts):
        """Tell our observes that we're about to change.
        """
        changeOpts['time'] = ChangeTimes.BEFORE
        self.__notifyList(self.__preObservers, key, **changeOpts)
        return
    
    def DidChangeValueForKey(self, key, **changeOpts):
        """Tell our observers that we have just changed.
        """
        changeOpts['time'] = ChangeTimes.AFTER
        self.__notifyList(self.__postObservers, key, **changeOpts)
        return
    
    def AddObserver(self, keyPath, callback, options=NotifyOptions()):
        """Make 'callback' an observer of changes to the attribute
        'attr'. The callback must be a callable taking one argument,
        which will be an Observer.Change object.
        """
        assert isinstance(keyPath, str)

        keyParts = keyPath.split('.', 1)

        # We always add an observer to the raw string
        self._AddObserverToLocalKey(keyPath, callback, options)
        if len(keyParts) == 2:
            localKey, restOfKey = keyParts
            self._AddObserverToComplexKey(localKey, restOfKey, callback, options)
            pass

        return
    
    def _AddObserverToComplexKey(self, localKey, restOfKey, callback, options):
        # Create the wrapper object
        wrapper = ComplexKeyCallbackWrapper(self, localKey, restOfKey, callback, options)
        # Make the wrapper observe the local part of the key path
        self._AddObserverToLocalKey(localKey, wrapper.LocalKeyPartChange,
                                    NotifyOptions(BEFORE_CHANGE=True,
                                                  AFTER_CHANGE=True))
        localValue = getattr(self, localKey)
        # And make it observe the sub-objects part of the key path
        localValue.AddObserver(restOfKey, wrapper.RestOfKeyChange, options)
        return
    
    def _AddObserverToLocalKey(self, keyPath, callback, options):
        for flag, obs in ((options.BEFORE_CHANGE, self.__preObservers),
                          (options.AFTER_CHANGE, self.__postObservers)):
            
            if flag:
                try:
                    attrObs = obs[keyPath]
                except KeyError:
                    attrObs = obs[keyPath] = set()
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
    
    def __getnewargs__(self):
        """Must ensure that __new__ is called before the state is set.
        Do this by returning an empty tuple of arguments.
        """
        return ()
    
    def __getstate__(self):
        picdic = {}
        for attr in self._Args:
            val = getattr(self, attr)
            if isinstance(val, ObservableList):
                picdic[attr] = [io for io in val]
            else:
                picdic[attr] = val
                pass
            continue
        return picdic
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        return
    
    def CloneFrom(self, other):
        try:
            attrList = copy(self._CloneOrder)
        except AttributeError:
            attrList = []
            pass
        
        for k in self._Args:
            if k not in attrList:
                attrList.append(k)
                pass
            continue
        
        for attr in attrList:
            val = getattr(self, attr)
            if isinstance(val, ObservableList):
                # first clear our list
                ourList = val
                while len(ourList):
                    ourList.pop()
                    continue
                # Copy in the new ones
                for obj in getattr(other, attr):
                    newObj = type(obj)()
                    newObj.CloneFrom(obj)
                    ourList.append(newObj)
                    continue
                
            elif isinstance(val, Observable):
                val.CloneFrom(getattr(other, attr))
            else:
                if(hasattr(other,attr)):
                    setattr(self, attr,
                            getattr(other, attr))
                pass
            continue
        return

    pass

class ComplexKeyCallbackWrapper(object):
    def __init__(self, obj, localKey, restOfKey, callback, options):
        self.obj = obj
        self.localKey = localKey
        self.restOfKey = restOfKey
        self.fullKey = localKey+'.'+restOfKey
        self.callback = callback
        self.options = options
        return

    def LocalKeyPartChange(self, change):
        # This is a change of the first part of the key
        if change.time is ChangeTimes.BEFORE:
            # Notify 
            self.obj.WillChangeValueForKey(self.fullKey)
            # Remove self from sub object's observers
            getattr(self.obj, self.localKey).RemoveObserver(self.restOfKey, self.RestOfKeyChange)
            pass
        if change.time is ChangeTimes.AFTER:
            # Notify 
            self.obj.DidChangeValueForKey(self.fullKey)
            # Add self as sub object's observer
            getattr(self.obj, self.localKey).AddObserver(self.restOfKey, self.RestOfKeyChange, self.options)
            pass
        return
    
    def RestOfKeyChange(self, change):
        # The rest of the key's value changed, here we just forward
        # the notification.
        if change.time is ChangeTimes.BEFORE:
            self.obj.WillChangeValueForKey(self.fullKey)
            pass
        if change.time is ChangeTimes.AFTER:
            self.obj.DidChangeValueForKey(self.fullKey)
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
            self.obj.WillChangeValueForKey(self.ofAttr)
            pass
        if change.time is ChangeTimes.AFTER:
            self.obj.DidChangeValueForKey(self.ofAttr)
            pass
        return
    pass

class ObservableList(Observable, collections.MutableSequence):
    """A list whose member mutations can be observed. Do this by
    observing special attributes '@INSERTION', '@REMOVAL',
    '@REPLACEMENT'.
    
    """
    _Args = {'_ObservableList__contents': []}
    def __init__(self, iterable=None):
        if iterable is None:
            self.__contents = list()
        else:
            self.__contents = list(iterable)
            pass

        self.AddDependency('length', '@INSERTION')
        self.AddDependency('length', '@REMOVAL')
        return
    
    def __copy__(self):
        return type(self)(self.__contents)
    
    def __len__(self):
        return self.__contents.__len__()
    
    def __getitem__(self, index):
        return self.__contents.__getitem__(index)
    
    def __setitem__(self, index, obj):
        if index >= len(self) or index < -len(self):
            raise IndexError('ObservableList assignment index out of range')
        
        self.WillChangeValueForKey('@REPLACEMENT', index=index)
        self.__contents.__setitem__(index, obj)
        self.DidChangeValueForKey('@REPLACEMENT', index=index)
    
    def __delitem__(self, index):
        if index >= len(self) or index < -len(self):
            raise IndexError('ObservableList assignment index out of range')
        
        self.WillChangeValueForKey('@REMOVAL', index=index)
        self.__contents.__delitem__(index)
        self.DidChangeValueForKey('@REMOVAL', index=index)
        return
    
    def insert(self, index, object):
        """Insert object before index.
        """
        self.WillChangeValueForKey('@INSERTION', index=index)
        self.__contents.insert(index, object)
        self.DidChangeValueForKey('@INSERTION', index=index)
        return

    @property
    def length(self):
        return self.__len__()
    
    def __str__(self):
        return self.__contents.__str__()
    def __repr__(self):
        return self.__contents.__repr__()
    pass

if __name__ == "__main__":
    import pdb
    
    class Observed(Observable):
        def __init__(self, x, pressure):
            self.x = x
            self.Pressure = pressure
            
            self.AddDependency('y', 'x')
            self.AddDependency('PressureString', 'Pressure.x')
            self.AddDependency('PressureString', 'Pressure.y')
            self.AddDependency('PressureString', 'Pressure.z')
            return
        
        @property
        def y(self):
            return self.x

        @property
        def PressureString(self):
            return 'p = %f + %f cos(wt + degtorad(%f))' % (self.Pressure.x,
                                                           self.Pressure.y,
                                                           self.Pressure.z)
        
        pass
    
    class Vector(Observable):
        def __init__(self, x,y,z):
            self.x = x
            self.y = y
            self.z = z
            return
        
    class Observer(object):
        def change(self, change):
            print change.key, change.obj.GetValueForKey(change.key)
            return
        pass

    print "Intantiate Observed a"
    a = Observed(7, Vector(1,2,3))
    b = Observer()
    print "Add b as observer of a.y, a.Pressure.x and a.PressureString"
    a.AddObserver('y', b.change)
    a.AddObserver('Pressure.x', b.change)
    a.AddObserver('PressureString', b.change)
    
    print "Set a.x"
    a.x = 1
    print "Set a.Pressure.x"
    a.Pressure.x = 3e8
    print "Set a.Pressure"
    a.Pressure = Vector(10,12,14)
    
    olist = ObservableList()
    print 'Add b as observer of olist.@INSERTION'
    olist.AddObserver('@INSERTION', b.change)
    print 'Append to olist'
    olist.append(6)
