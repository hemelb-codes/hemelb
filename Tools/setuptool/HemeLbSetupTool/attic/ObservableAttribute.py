import pdb
import collections
from Enum import Enum

class NoSuch(object):
    """Special None-like object for when an attribute did not exist.
    """
    pass
NoSuch = NoSuch()
ChangeTimes = Enum("BEFORE", "AFTER")

class NotifyOptions(object):
    def __init__(self, **kwargs):
        defaults = {'BEFORE_CHANGE': False,
                    'AFTER_CHANGE': True}
        
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
        if time not in ChangeTimes:
            raise ValueError("Keyword argument 'time' must be specified and one of %s" %
                             str(ChangeTimes))
        self.time = time
        self.obj = obj
        self.new = new
        self.old = old
        self.index = index
        return

    pass

class HasObservableProperties(object):
    """Base class for objects that will need to be observed.
    """
    
    def __new__(cls, *args, **kwargs):
        # Create the instance as normal
        new = object.__new__(cls)

        # Check class dict for Observables
        for name, val in cls.__dict__.iteritems():
            if isinstance(val, Observable):
                # Init the observable for this instance
                val._Init(name, new, cls)
                pass
            continue
        for name, val in cls.__dict__.iteritems():
            if isinstance(val, Observable):
                # Init the observable for this instance
                val._AddDependencyObservers(new)
                pass
            continue
        return new
    
    def AddObserver(self, propertyName, callback, options=NotifyOptions()):
        """Make 'callback' an observer of changes to propertyName. The
        callback must be a callable taking one argument, which will be
        an Observer.Change object.
        """
        return getattr(self.__class__, propertyName).AddObserver(self, callback, options)
    
    def RemoveObserver(self, propertyName, callback):
        """Remove an observer.
        """
        return getattr(self.__class__, propertyName).RemoveObserver(self, callback)
    
    def AddDependency(self, depending, controlling, instance=None):
        if instance is None:
            instance = self
            pass
        
        def closure(change):
            if change.time is ChangeTimes.BEFORE:
                getattr(self.__class__, depending).WillChange(instance)
                pass
            if change.time is ChangeTimes.AFTER:
                getattr(self.__class__, depending).DidChange(instance)
                pass
            return
        
        self.AddObserver(depending, closure,
                         options=NotifyOptions(BEFORE_CHANGE=True,
                                               AFTER_CHANGE=True))
        return

    def WillChangeAttr(self, attr, **opts):
        self.__class__.__dict__[attr].WillChange(self, **opts)
        return
    def DidChangeAttr(self, attr, **opts):
        self.__class__.__dict__[attr].DidChange(self, **opts)
        return
    
    #     if instance is None:
    #         instance = self
    #         pass
    #     return getattr(self.__class__, depending).AddDependency(self,
    #                                                             controlling,
    #                                                             instance)
    pass

class Observable(object):
    def __init__(self, *args, **kwargs):
        self.dependencies = set()
        return
    
    def _Init(self, name, instance, cls):
        """Private method. To be called from the containing class
        (i.e. HasObservableProperties or subclass thereof) during
        __new__.

        Sets up anything necessary on the instance, mainly the shadow
        dictionary use to store our private data.
        """
        self.name = name
        self.shadowName = '__%s_%s_shadow' % (cls.__name__,
                                              name)
        sdict = {}
        sdict['preObservers'] = set()
        sdict['postObservers'] = set()
        
        setattr(instance, self.shadowName, sdict)
        
        return
    
    def _AddDependencyObservers(self, instance):
        """Private method. Add the observers necessary for dependency
        handling.
        """
        for prop in self.dependencies:
            def closure(change):
                if change.time is ChangeTimes.BEFORE:
                    self.WillChange(instance)
                    pass
                if change.time is ChangeTimes.AFTER:
                    self.DidChange(instance)
                    pass
                return
            
            prop.AddObserver(instance, closure,
                             options=NotifyOptions(BEFORE_CHANGE=True,
                                                   AFTER_CHANGE=True))
            continue
        return
    
    def GetShadowDict(self, instance):
        """Private method. Get the shadow dict from the instance.
        """
        return getattr(instance, self.shadowName)

    def _NotifyList(self, instance, oSet, attr, **changeOpts):
        opts = dict(obj=instance,
                    attr=attr)
        opts.update(changeOpts)
        # Construct a change object
        change = Change(**opts)
        
        for o in oSet:
            # Call each observer callback with the change object
            o(change)
            continue
        return

    def WillChange(self, instance, **changeOpts):
        """Tell our observers on instance that we're about to change.
        """
        changeOpts['time'] = ChangeTimes.BEFORE
        self._NotifyList(instance,
                          self.GetShadowDict(instance)['preObservers'],
                          self.name, **changeOpts)
        return
    
    def DidChange(self, instance, **changeOpts):
        """Tell our observers on instance that we have just changed.
        """
        changeOpts['time'] = ChangeTimes.AFTER
        self._NotifyList(instance,
                          self.GetShadowDict(instance)['postObservers'],
                          self.name, **changeOpts)
        return
    
    def AddObserver(self, instance, callback, options=NotifyOptions()):
        """Make 'callback' an observer of our changes on instance. The
        callback must be a callable taking one argument, which will be
        an Observer.Change object.
        """
        
        for flag, obs in ((options.BEFORE_CHANGE,
                           self.GetShadowDict(instance)['preObservers']),
                          (options.AFTER_CHANGE,
                           self.GetShadowDict(instance)['postObservers'])):
            
            if flag:
                obs.add(callback)
                pass
            continue
        
        return

    def RemoveObserver(self, instance, observer):
        """Remove the observer from the list to be notified on change.
        """
        for timeObservers in (self.GetShadowDict(instance)['preObservers'],
                              self.GetShadowDict(instance)['postObservers']):
            try:
                timeObservers.remove(observer)
            except KeyError:
                pass
            continue
        
        return

    def AddDependency(self, onProperty):
        """Mark this property as dependent on the value of another
        (onProperty). Thus when onProperty is changed, the observers
        of self will be notified.
        """
        
        if onProperty in self.dependencies:
            # Trying to add a redundant dependency
            return
        self.dependencies.add(onProperty)
        
        return

    pass

class ObservableProperty(Observable, property):
    """ObservableProperty(fget=None, fset=None, fdel=None, doc=None)-> property attribute
    
    fget is a function to be used for getting an attribute value, and
    likewise fset is a function for setting, and fdel a function for
    del'ing, an attribute.  Typical use is to define a managed
    attribute x:
    
    class C(object):
        def getx(self): return self._x
        def setx(self, value): self._x = value
        def delx(self): del self._x
        x = property(getx, setx, delx, \"I'm the 'x' property.\")
    
    Decorators make defining new properties or modifying existing ones
    easy:

        class C(object):
            @property
            def x(self): return self._x
            @x.setter
            def x(self, value): self._x = value
            @x.deleter
            def x(self): del self._x
    """
    
    def __init__(self, *args, **kwargs):
        """See class __doc__.
        """
        Observable.__init__(self, *args, **kwargs)
        property.__init__(self, *args, **kwargs)
        return
    
    def __set__(self, instance, value):
        """Implementation of data descriptor mechanism. See
        http://docs.python.org/reference/datamodel.html#descriptors
        Wraps the property method, triggering pre- and post-
        notifications.
        """
        if self.fset is not None:
            self.WillChange(instance)
            
        property.__set__(self, instance, value)
        self.DidChange(instance)
        return
        
    pass

class ObservableAttribute(ObservableProperty):
    """Makes a property that acts just like a normal attribute, but
    notifies its observers of changes.
    """
    def __init__(self, default=NoSuch):
        self.default = default
        ObservableProperty.__init__(self, self.getattr, self.setattr, self.delattr)
        pass
    
    def _Init(self, name, instance, cls):
        ObservableProperty._Init(self, name, instance, cls)
        self.GetShadowDict(instance)['value'] = self.default
        return
    
    def getattr(self, instance):
        try:
            ans = self.GetShadowDict(instance)['value']
            if ans is NoSuch: raise KeyError
            
        except KeyError:
            raise AttributeError("'%s' object has no attribute '%s'" %
                                 (ownerClass.__name__, self.name))
        return ans
    
    def setattr(self, instance, value):
        self.GetShadowDict(instance)['value'] = value
        return
    
    def delattr(self, instance):
        try:
            del self.GetShadowDict(instance)['value']
        except KeyError:
            raise AttributeError(self.name)
        return
    pass

class ObservableAction(Observable):
    def __init__(self, action=None):
        self.action = None
    pass

class ObservableList(HasObservableProperties, collections.MutableSequence):
    """A list whose member mutations can be observed. Do this by
    observing special attributes '@INSERTION', '@REMOVAL',
    '@REPLACEMENT'.
    
    """
    pdb.set_trace()
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

setattr(ObservableList, '@REPLACEMENT', ObservableAction())
setattr(ObservableList, '@REMOVAL', ObservableAction())
setattr(ObservableList, '@INSERTION', ObservableAction())
        
if __name__ == "__main__":
    import pdb
    class O(HasObservableProperties):
        x = ObservableAttribute()
        y = ObservableProperty(lambda self: self.x)
        y.AddDependency(x)
        
        pass
    
    def cb(change):
        print change.__dict__
        return
    
    o = O()
    o.AddObserver('y', cb)
    o.x = 5
    
    olist = ObservableList()
    olist.AddObserver('@INSERTION', cb)
    olist.append(4)
