# HemeLB

A typical workflow with HemeLB consists of four steps:

1. A preprocessing step where you create a mesh or geometry file (we
   use extension .gmy) from a representation of the surface of your
   domain. For this, you need to use the application in the
   [geometry-tool](../geometry-tool)
   directory. [Documentation](user/geometry-tool.md)

2. Setting options in the configuration XML file, such as timestep and
   output data and
   frequency. [Documentation](users/XmlConfiguration.md)

3. Compiling and running the
   [main application](user/main-application.md).

4. Post processing the results, using the Python packages in the
   [python-tools](../python-tools) directory. Documentation TBC


[Developer documentation](dev)

