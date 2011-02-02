from wx.lib.pubsub import Publisher

class Notifier(object):
    """This base class will send a pubsub message to a specified topic
    when a specified attribute has its value changed.
    """
    
    def __new__(cls, *args, **kwargs):
        """Do required initialisation to ensure that we insert the
        private attribute that stores the topics that will be messaged
        on change.
        """
        new = object.__new__(cls)
        object.__setattr__(new, '_topicsToNotify', {})
        return new
    
    @classmethod
    def ClassTopic(cls, subtopic=None):
        """Get a dot-delimited string that describes the instance's
        class. Usable as a topic string.
        """
        topic = cls.__module__+ '.' + cls.__name__
        if subtopic is not None:
            topic += '.' + subtopic
        return topic
        
    def SelfTopic(self, subtopic=None):
        """Create a dot-delimited string that uniquely (within the
        process) indentifies a particular instance.
        """
        topic = self.ClassTopic(str(id(self)))
        if subtopic is not None:
            topic += '.' + subtopic
        return topic
        
    def __setattr__(self, name, value):
        """Set the attribute and notify topics if the value's changed.
        """
        try:
            # If the value hasn't changed, don't need to notify 
            if value == getattr(self, name):
                return
        except AttributeError:
            # But if this is the first assignment to the attribute,
            # the above will be an error; one that we can swallow.
            pass
        
        try:
            tList = self._topicsToNotify[name]
        except KeyError:
            # No topics to notify, set attr and return
            return object.__setattr__(self, name, value)
        
        for topic in tList:
            # Send the notifications to the topics
            Publisher.sendMessage(topic, data=value)
            continue
        # Set attribute and return
        return object.__setattr__(self, name, value)
    
    def addNotifiedTopic(self, attr, topic):
        """Append a topic to the list to be notified when the
        attribute changes. Returns a boolean indicating whether the
        topic was added or not.
        """
        assert isinstance(attr, str)
        assert hasattr(self, attr)
        
        topic = self._canonicaliseTopic(topic)
        try:
            attrTopics = self._topicsToNotify[attr]
        except KeyError:
            attrTopics = self._topicsToNotify[attr] = []
            pass
        
        if topic in attrTopics:
            return False
        else:
            attrTopics.append(topic)
            return True
        return

    def removeNotifiedTopic(self, attr, topic):
        """Remove a topic from the list to be notified on
        change. Returns a boolean indicating whether the topic was
        removed or if it wasn't present.
        """
        assert isinstance(attr, str)
        assert hasattr(self, attr)
        topic = self._canonicaliseTopic(topic)
        
        try:
            attrTopics = self._topicsToNotify[attr]
        except KeyError:
            return False

        try:
            attrTopics.remove(topic)
            return True
        except ValueError:
            return False
        return
    
    @staticmethod
    def _canonicaliseTopic(topic):
        if isinstance(topic, str):
            return topic
        else:
            return '.'.join(topic)
        return
    
    pass

if __name__ == "__main__":
    class Observed(Notifier):
        def __init__(self):
            self.x = 0
            return
        pass
    
    class Observer(object):
        def changex(self, msg):
            print msg.topic, msg.data
            return
        pass
    
    a = Observed()
    
    b = Observer()
    
    a.addNotifiedTopic('x', 'changex')
    Publisher.subscribe(b.changex, 'changex')
    
    a.x = 1
    
