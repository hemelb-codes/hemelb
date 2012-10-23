# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

import os.path
import cPickle
from functools import wraps

def cache(filename, *argPatterns, **kwargPatterns):
    """Given a name of a file, returns a decorator which will cache
    the results of a function call (of a function with a first arg
    of "base") to "filename"
    """
    def decorator(f):
        # this is the actual decorator
        # f is the function to be modified
        
        @wraps(f)
        def cacher(*args, **kwargs):
            # create the full path to the cache file
            cachefile = filename
            if len(argPatterns):
                cachefile += '.'+'.'.join([a%args[i] 
                                           for i,a in enumerate(argPatterns)])
                pass
            
            keys = kwargPatterns.keys()
            if len(keys):
                keys.sort()
                cachefile += '.'+'.'.join([kwargPatterns[k]%kwargs[k] 
                                           for k in keys])
                pass
            
            
            if os.path.exists(cachefile):
                ans = cPickle.load(file(cachefile))
            else:
                ans = f(*args, **kwargs)
                cPickle.dump(ans,
                             file(cachefile, 'wb'),
                             protocol=2)
                pass
            return ans
        
        # return the wrapped function
        return cacher

    return decorator

def processesFile(f):
    """Decorator which will cache the results of a function call (of a
    function with a single argument that specifies a file which is the
    sole input) to "outfile".
    """
    
    @wraps(f)
    def cacher(infile):
        # create the full path to the cache file
        cachefile = infile + '.' + f.__module__ + '.' + f.func_name + '.cache'
        
        if os.path.exists(cachefile) and \
                os.path.getmtime(cachefile) > os.path.getmtime(infile):
            # Is the cachefile newer than the input?
            ans = cPickle.load(file(cachefile))
        else:
            ans = f(infile)
            cPickle.dump(ans,
                         file(cachefile, 'wb'),
                         protocol=2)
            pass
        return ans
    
    # return the wrapped function
    return cacher

