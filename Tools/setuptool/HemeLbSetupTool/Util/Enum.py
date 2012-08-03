# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

class Enum(object):
    def __init__(self, *args):
        """Create an Enum with the list of names given. Each item must
        either be a string or a pair whose first value is a string and
        second is an int greater than all preceeding ints.
        """
        
        class EnumItem(object):
            def __init__(slf, i, name):
                assert name.isalnum() and name[0].isalpha()
                slf.i = i
                slf.name = name
                slf.enum = self
                return

            def __repr__(slf):
                return "<EnumItem " + slf.name + " = %d>" % slf.i
            pass
        
        self.ItemClass = EnumItem
        
        i = 0
        names = []
        values = []
        for item in args:
            if isinstance(item, str):
                names.append(item)
                values.append(i)
            else:
                name, i = item
                assert isinstance(name, str)
                assert isinstance(i, int)
                assert len(values)==0 or values[-1] < i
                
                names.append(name)
                values.append(i)
                pass
            i += 1
            continue
        
        self.itemsByIndex = items = tuple(EnumItem(*arg) for arg in zip(values, names))
        self.itemsByName = dict((item.name, item) for item in items)
        
        return
    
    def __getitem__(self, i):
        return object.__getattribute__(self, 'itemsByIndex')[i]
    
    def __getattribute__(self, attr):
        try:
            return object.__getattribute__(self, 'itemsByName')[attr]
        except KeyError:
            raise AttributeError("No such constant '%s' in enum" % attr)
        
    def __contains__(self, item):
        items = object.__getattribute__(self, 'itemsByIndex')
        return items.__contains__(item)

    def __repr__(self):
        ans = "<Enum "
        first = True
        for item in object.__getattribute__(self, "itemsByIndex"):
            if not first: ans += ","
            ans += " " + item.name + " = %d" % item.i
            first = False
            continue
        ans += ">"
        return ans
    
    pass
