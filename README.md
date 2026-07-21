# TwoPunctures (GRChombo/GRTeclyn)

[![GRChombo](https://github.com/GRTLCollaboration/TwoPunctures/actions/workflows/grchombo.yml/badge.svg)](https://github.com/GRTLCollaboration/TwoPunctures/actions/workflows/grchombo.yml)
[![GRTeclyn](https://github.com/GRTLCollaboration/TwoPunctures/actions/workflows/grteclyn.yml/badge.svg)](https://github.com/GRTLCollaboration/TwoPunctures/actions/workflows/grteclyn.yml)
![Style](https://github.com/GRChombo/TwoPunctures/workflows/Style/badge.svg)

This repository originates from the
[standalone TwoPunctures](https://bitbucket.org/relastro/twopunctures-standalone)
code which in turn comes from the original Einstein Toolkit thorn
[TwoPunctures](https://bitbucket.org/einsteintoolkit/einsteininitialdata/).

It has been modified to allow integration with
[GRChombo](https://github.com/GRTLCollaboration/GRChombo) and
[GRTeclyn](https://github.com/GRTLCollaboration/GRTeclyn).

##  Physics

The code creates Bowen-York initial data for two puncture black holes using a
single domain spectral method. The method is described in _Marcus Ansorg,
Bernd Brügmann, Wolfgang Tichy_, "A single-domain spectral method for black hole
puncture data"
[PRD 70, 064011 (2004)](https://doi.org/10.1103/PhysRevD.70.064011)
([arXiv:gr-qc/0404056](https://arxiv.org/abs/gr-qc/0404056)).

## GRChombo integration

To build the BinaryBH GRChombo example with TwoPunctures initial data, set the
environment variable `TWOPUNCTURES_SOURCE` to the location of the Source
subdirectory of this repository. Note that the parameters used are different to
the vanilla BinaryBH example.

There is also a standalone example in this repository which can be used to check
the bare and ADM masses of the final TwoPunctures result. To build this code,
the `GRCHOMBO_SOURCE` environment variable must be set to the Source
subdirectory of the GRChombo repository.

In either case linking with the GNU Scientific Library (GSL) is required.

## GRTeclyn integration

To build the BinaryBH GRTeclyn example with TwoPunctures initial data, set the
makefile variable `USE_TWOPUNCTURES = TRUE`. By default, the build system will
assume this repository is cloned next to GRTeclyn i.e.
`/path/to/GRTeclyn/../TwoPunctures` (like AMReX) but this can be overridden by
setting the `TWOPUNCTURES_HOME` variable. Both of these variables can be set in
the following ways:

* on the command line e.g.
  ```
  make -j 4 USE_TWOPUNCTURES=TRUE TWOPUNCTURES_HOME=/path/to/TwoPunctures
  ```
* in the `Make.local-pre` or `Make.local` files located in `amrex/Tools/GNUMake`
* by setting it as an environment variable

As for GRChombo, it requires linking with the GNU Scientific Library (GSL). If
necessary, the GSL install location can be passed to the build system by setting
`GSL_HOME`.

> [!NOTE]
> TwoPunctures parameter names differ between GRChombo and GRTeclyn.

> [!NOTE]
> When running GRTeclyn on GPUs, TwoPunctures will still run and
> interpolate the data onto the grid on the CPU. In this case, since you will
> generally have a much greater number of cells per MPI process, spectral
> interpolation will be *very* slow.

## License

The authors of the original TwoPunctures code are Marcus Ansorg, Erik Schnetter,
Frank Löffler. The modification to the standalone code on which this repo is
based was done by Federico Guercilena, Sven Köppel

This code is licensed under the LGPLv2.1 license. Please see the LICENSE file
for details.
