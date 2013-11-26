class default_property(object):
    """Decorator to make a property with a default value (the result
    of evaluating the method) that can be overriden by setting a
    value. The default is restored by deleting the attribute.
    """
    def __init__(self, fcalc):
        self.value_attr_name = '_default_' + fcalc.func_name
        self.calc_default = fcalc
        return
    
    def __get__(self, instance, owner_cls):
        try:
            return getattr(instance, self.value_attr_name)
        except AttributeError:
            return self.calc_default(instance)

    def __set__(self, instance, value):
        return setattr(instance, self.value_attr_name, value)

    def __delete__(self, instance):
        try:
            return delattr(instance, self.value_attr_name)
        except AttributeError:
            raise AttributeError("can't delete attribute")
        
    pass
