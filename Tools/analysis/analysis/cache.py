import os.path
import cPickle

class cache_property(object):
    """Decorator for caching computed properties to disc.  Only
    defined for classes which are duck-typed like
    HemeLbRunResults. Results are cached in the run directory as
    $RunDirectory/$PropertyName.cache with the pickle format.
    """
    def __init__(self, fget):
        self.getter = fget
        self.cache_file_basename = fget.func_name + '.cache'

    def __get__(self, instance, owner_cls):
        cache_file = os.path.join(instance.RunDirectory, self.cache_file_basename)

        if os.path.exists(cache_file):
            with file(cache_file) as cf:
                return cPickle.load(cf)

        ans = self.getter(instance)
        with file(cache_file, 'w') as cf:
            cPickle.dump(ans, cf, protocol=cPickle.HIGHEST_PROTOCOL)
        return ans
