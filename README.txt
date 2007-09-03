HemeLB Version 1.3-RG

0. Overview
1. Compiling HemeLB
2. Running HemeLB
3. RealityGrid Steering and Visualisation
4. Ray Tracing Outline




0. Overview ==============================================================

HemeLB's documentation is in "doc/HemeLB_doc.pdf". It comprises the paper
submitted to Comput. Phys. Comm. (with the same title and authors) and
the usage at the end of the file. The code is in directory hemelb/Code, two
input configurations and parameters are in the directory
hemelb/Input" (see the documentation for details).

Please feel free to email me at m.mazzeo@ucl.ac.uk if you have any
problems.

Please email the GENIUS discussion mailing list, genius-discuss@ucl.ac.uk
if you have any problems, or wish to report bugs, etc.



1. Compiling HemeLB ======================================================

HemeLB can be run in standalone mode or RealityGrid steering mode.
The definition RG is used to seperate out the RealityGrid specific parts
of HemeLB.



2. Running HemeLB ========================================================

HemeLB is run with a single argument that gives the location of an input
file.

Two sets of input/output files can be found in hemelb/Input. The current
version of hemelb requires full path names to be specified in these files

Input/square_duct_32x16x16_input.asc
Input/angio1_input.asc

Before running HemeLB please edit these files appropriately.


3. RealityGrid Steering and Visualisation =============================


4. Ray Tracing outline ================================================

The flow field is rendered by a parallelised raytracer which is built into
the same code as the LB-fluid solver. In other words, the rendering of the 
HemeLB is not done externally, which is usually the case, but within the same
application, exploiting the domain decomposition of the LB solver. The code 
to perform the ray tracing can be found in Code/rt.cc. Technical details on
the ray tracing core and the parallel approach implemented are provided
elsewhere.

The parameters used by the ray-tracing algorithm are read from a file
whose name (the path should be specified if not defined in a RSL
script) must occupy the 4-th line of the input file which is read on
command line as first argument when submitting a job (see
hemelb/Doc/HemeLB_doc.pdf). The rendering parameters must be inserted in
different lines and are

pixels_x
pixels_y
ctr_x
ctr_y
ctr_z
longitude
latitude
zoom
image_frequency
flow_field_type
is_isosurface
absorption_factor
cutoff
density_max
velocity_max
stress_max

"pixels_x" and "pixels_y" are the pixels along the X axis and Y one
respectively; "ctr_x", "ctr_y", "ctr_z" specify the center of the
scene. The viewpoint looks towards that point from the viewpoint
angles "longitude" and "latitude" which define the polar and vertical
angles respectively. "zoom" changes the projection size of the
rendering.  The latter is performed at every "image_frequency" LB time
steps. "flow_field_type" must be 0, 1, or 2 if the pressure/density,
velocity or effective von Mises stress flow fields respectively must
be rendered. The iso-surface is calculated if "is_isosurface" is 1,
otherwise the volume rendering with "absorption_factor" as opacity
factor is performed. The parameter "cutoff" must be between 0 and 1;
if the velocity field is on, the volume rendering (or the iso-surface)
is calculated considering only the fluid lattice sites which have a
velocity greater than or equal to cutoff * velocity_max. A similar
argument holds for the density/pressure and effective stress rendering
modes.

The image is outputted on file which must be specified in the 5-th line
of the input file which is read on command line as first argument when
submitting a job. Its format is under investigation and will be
descripted elsewhere.



This document was last updated on the 3/9/2007.

