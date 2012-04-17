#Jiggery-pokery so that client classes of RemoteHemeLB can write remote.Latitude=50

class SteeredParameter(object):
    def __init__(self,index,name):
        self.name=name
        self.index=index
    def initialise_in_instance(self,instance,value):
        if not hasattr(instance,'properties'):
            instance.properties={}
            instance.property_old_values={}
        if not self.name in instance.properties:
            instance.properties[self.name]=value
            instance.property_old_values[self.name]=value
    def changed(self,instance):
        return not self.old_value(instance)==self.value(instance)
    def old_value(self,instance):
        return instance.property_old_values[self.name]
    def value(self,instance):
        return instance.properties[self.name]
    def __get__(self,instance,owner):
        return self.value(instance)
    def __set__(self,instance,value):
        instance.property_old_values[self.name]=instance.properties[self.name]
        instance.properties[self.name]=value

    