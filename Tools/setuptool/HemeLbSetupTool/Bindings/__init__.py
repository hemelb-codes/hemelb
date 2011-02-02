"""This module contains a simple implementation of the Bindings
pattern for MVC applications, using that in Cocoa as its
inspiration. See:

http://developer.apple.com/library/mac/documentation/Cocoa/Conceptual/CocoaBindings/CocoaBindings.html

for details.

Basically your Model object(s) must inherit from Observable and have
public attributes/properties (my convention is to start these with an
upper case letter) that will be bound to View items. It can also contain methods that we

Your Controller should be an ObjectController or a subclass.


"""
import Bindings
import Mappers
import Translators
import Controllers

if __name__ == "__main__":
    from test import test
    test()
    
