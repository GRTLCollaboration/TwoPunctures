# TwoPunctures (GRChombo)

![Build BinaryBH example with TwoPunctures](https://github.com/GRChombo/TwoPunctures/workflows/Build%20BinaryBH%20example%20with%20TwoPunctures/badge.svg)
![Check Clang Format](https://github.com/GRChombo/TwoPunctures/workflows/Check%20Clang%20Format/badge.svg)

This repository originates from the
[standalone TwoPunctures](https://bitbucket.org/relastro/twopunctures-standalone)
code which in turn comes from the original Einstein Toolkit thorn
[TwoPunctures](https://bitbucket.org/einsteintoolkit/einsteininitialdata/).

It has been modified to allow integration with GRChombo.

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

## License

The authors of the original TwoPunctures code are Marcus Ansorg, Erik Schnetter,
Frank Löffler. The modification to the standalone code on which this repo is
based was done by Federico Guercilena, Sven Köppel

This code is licensed under the LGPLv2.1 license. Please see the LICENSE file
for details.
