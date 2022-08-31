# HemeLB

A typical workflow with HemeLB consists of four steps:

1. A preprocessing step where you create a mesh or geometry file (we
   use extension .gmy) from a representation of the surface of your
   domain. For this, you need to use the application in the
   [geometry-tool](../geometry-tool)
   directory. [Documentation](user/geometry-tool.md)

2. Setting options in the configuration XML file, such as timestep and
   output data and
   frequency. [Documentation](user/XmlConfiguration.md)

3. Compiling and running the
   [main application](user/main-application.md).

4. Post processing the results, using the Python packages in the
   [python-tools](../python-tools) directory. [Documentation](user/python-tools.md)


[Developer documentation](dev)

# Publications

## Key code publications

- R.W. Nash, H.B. Carver, M.O. Bernabeu, J. Hetherington, D. Groen, T.
  Krüger, P.V. Coveney, "Choice of boundary condition for
  lattice-Boltzmann simulation of moderate-Reynolds-number flow in
  complex domains", Phys. Rev. E (2014).
  https://doi.org/10.1103/PhysRevE.89.023303

- D. Groen, J. Hetherington, H.B. Carver, R.W. Nash, M.O. Bernabeu,
  "Analysing and modelling the performance of the HemeLB
  lattice-Boltzmann simulation environment", J. Comput. Sci. (2013).
  https://doi.org/10.1016/j.jocs.2013.03.002

- M.D. Mazzeo & P.V. Coveney, "HemeLB: A high performance parallel
  lattice-Boltzmann code for large scale fluid flow in complex
  geometries", Comput. Phys. Commun. (2008)
  https://doi.org/10.1016/j.cpc.2008.02.013


## Papers using HemeLB
