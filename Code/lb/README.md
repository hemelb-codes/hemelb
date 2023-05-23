# Description of the LB namespace

In `concepts.h` we give C++20 concepts specifying exactly what members etc.
need to be present. Here we explain at a higher level.

## Lattice (or velocity set)

E.g. D3Q15 etc. Describes the set of lattice velocity vectors used for the LB
and provides functions to calculate hydrodynamic variables from the
distribution functions, and equilibrium distribution functions from those.

These are implemented as a class template - see lattices/Lattice.h.

## Kernel

These implement most of the collision part of the LB algorithm. Requires a
lattice as a template parameter. Examples are LBGK, MRT, entropic, etc.
Does operations such as calculate the hydrodynamic variables and the post-
collisional distributions.

Importantly, this type must have a public member type (`VarsType`) that
provides/refers to storage for the variables needed such as the initial
distributions, the relevant hydrodynamic variables, and the post-collision
distributions. This should be an instantiation of `HydroVars<KernelType>`,
and may be specialised if extra fields are needed (e.g. MRT).

## Collision

This gives the collision the opportunity to override the behaviour of the
collision kernel. This was used in the original Mazzeo & Coveney papers
to impose wall and inlet/outlet boundary conditions. Currently, only the
`Normal` collision is used in the main code.

## Streamer

These types implement the combined collide-and-stream operation over a
contiguous range of sites. The design is to for the streamers to be pretty
simple and for them to basically just control iteration over the sites and
directions while delegating the logic of actually streaming to other classes
along the links (i.e. the velocity set). This works well for all boundary
conditions except Junk-Yang.

We distinguish three types of link streamers:

- bulk
- wall
- iolet
there being only one bulk link type but multiple choices for wall and iolet.
These are combined by the StreamerTypeFactory template (with a `NullLink`
type used to indicate the absence of that type) to produce a streamer.

Junk-Yang cannot work per-link so has to be special cased.

## Headers

The headers `Collisions.h`, `Kernels.h`, `Lattices.h`, `Streamers.h` provide
access to all the types fulfilling the respective concept, as well as an
alias for the default collision etc that has been configured by CMake.

## LbmParameters and InitParams

These types define the relaxation towards equilibrium etc and give access
into all the parameters used to initialise a simulation.

## LBM<Traits>

This class combines the various streamers (via the `Traits` template class) and
runs them in the right order for the communications, i.e.:
1. post receives,
2. do domain edge sites,
3. asynchronously send their data,
4. do mid-domain sites,
5. wait on communications,
6. fix up edge sites

