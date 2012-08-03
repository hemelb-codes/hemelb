# 
# Copyright (C) University College London, 2007-2012, all rights reserved.
# 
# This file is part of HemeLB and is CONFIDENTIAL. You may not work 
# with, install, use, duplicate, modify, redistribute or share this
# file, or any part thereof, other than as allowed by any agreement
# specifically made by you with University College London.
# 

from . import Translators

class ValueBinding(object):
    """This object manages the binding of zero or more controls to a
    model variable. It intermediates between the model mapper and
    widget mapper(s).
    """
    def __init__(self, modelMapper):
        modelMapper.SetBinding(self)
                
        self.modelMapper = modelMapper
        self.widgetMappers = set()
        self.updateSource = None
        
        modelMapper.Observe()
        
        return

    def GetAllMappers(self):
        return self.widgetMappers.union((self.modelMapper,))
    
    def BindWidget(self, widgetMapper):
        widgetMapper.SetBinding(self)
        
        self.widgetMappers.add(widgetMapper)
        self.Sync(self.modelMapper, widgetMapper)
        widgetMapper.Observe()
        return
    
    def MapperWasUpdated(self, source):
        """One of our mappers got triggered.
        """
        # Check we're not already responding to an update
        if self.updateSource is not None:
            return
        
        try:
            # Note who started this update
            self.updateSource = source
            
            # We're going to construct an iterable of the arguments to
            # Sync. Handle two cases differenly
            if source is self.modelMapper:
                # Just syncs from the model to all widgets
                syncs = ((source, wm) for wm in self.widgetMappers)
            elif source in self.widgetMappers:
                # Want to first sync widget->model, then model->other
                # widgets. First put all the syncs in the queue
                syncs = [(self.modelMapper, wm)
                         for wm in self.widgetMappers
                         if wm is not source]
                # Then have widget->model jump in front
                syncs.insert(0, (source, self.modelMapper))
            else:
                # it's not one of the mappers for this instance
                raise ValueError('Mapper is not one of mine!')

            # execute the syncs
            for args in syncs:
                self.Sync(*args)
                continue
            
        finally:
            # Make sure we note that the update's finished.
            self.updateSource = None
            pass
        return
    
    def Sync(self, source, dest):
        try:
            dest.Set(source.Get())
        except Translators.FormattingError:
            # Can't format this properly, so just give up for now
            pass
        except Translators.ValidationError:
            # We got a value from the widget (i.e. source) that
            # doesn't convert to the model requirements. Put back the
            # last known good
            source.Set(self.modelMapper.Get())
            pass
        
        return
    
    pass

class ActionBinding(object):
    pass

class WxActionBinding(ActionBinding):
    def __init__(self, widget, event):
        self.widget = widget
        self.event = event
        return

    def Bind(self, callback):
        self.widget.Bind(self.event, callback)
        return
    pass
