# GWHL
repository for the GWHL

The "Gravitational Wave Helper Lib" is a very basic math library with flat data types and code as simple as possible.

Some basic numeric routines are present:
- basic 3d vector BLAS 1,2,3, with non-standard interfaces
- a basic Cholesky solver (described in Lindegren, L.; Lammers, U.; Hobbs, D.; O'Mullane, W.; Bastian, U.; Hernández, J.: "The astrometric core solution for the Gaia mission", Appendix C [https://arxiv.org/abs/1112.4139])
- extremely fast and stable computation of associated Legendre polynomials 
- unit conversion functions
- some quaternion routines (rotation, normalization, ...)
- routines to compute spherical harmonics
- rotation matrices
- routines computing vector spherical harmonics (see F. Mignard and S. Klioner: Analysis of astrometric catalogues with vector spherical harmonics [https://www.aanda.org/articles/aa/abs/2012/11/aa19927-12/aa19927-12.html])

All routines have in common that they work on most basic data-types: linear double arrays instead of a special struct for vectors or matrices. That makes the code quite "compiler friendly" and hence fast.

This library is work in progress: no API stability is guaranteed, features can be added and removed at any time, no completeness nor correctness is guaranteed.

Acknowledgements / Contributors:
  Sergei A. Klioner (GWHL implements his mathematical definition of VSH routines)
  Lennart Lindegren (GWHL implements his modified cholesky solver)
  
  This library has been partly developed during projects financially supported by
ESA under grant 4000115263/15/NL/IB and the BMWi grant 50QG1402 awarded by the DLR.

## installation
* Get a copy of GNU GSL ([https://www.gnu.org/software/gsl/]), either by manual installation or via your distribution.
* Set the GSLDIR variable in the GWHL Makefile to this installation of GSL
* Adapt the rest of the Makefile to your needs
* make
* the resulting libGWHelper.so can be used as needed
