Documentation of glob_stab               {#mainpage}
==========================
This is an implementation of the stability analysis algorithm presented in
\cite Schauss2017. The method is explained in the following section.

Abstract from paper
-------------------
A novel approach to stability analysis of linear time-invariant (LTI) time-delay
systems of retarded type with incommensurate time delays and polynomial
dependence on parametric uncertainties is presented. Using a branch and bound
algorithm, which relies on Taylor Models and polynomials in Bernstein form, we
first determine the stability-crossing set in the delay/parameter space and then
evaluate stability for each disjoint region. The novel approach can be used to
either non-conservatively check stability of time-delay systems with interval
parameters or to map stable regions to a low-dimensional parameter space while
taking additional interval parameters into account. The algorithm is applied to
several examples with parametric uncertainties and/or incommensurate time
delays.

Software Structure
------------------
This implementation aims at providing reusable modules, especially for the
low-level mathematics. The overall software is split into four libraries:

* math:
  Consists of low-level math, especially different data types, e.g.,
    * Interval arithmetics (Iv)
    * Multivariate polynomials (Polynomial)
    * Taylor Models (TaylorModel)
    * Bernstein polynomials
    * Complex Bernstein polynomials (Bernstein2d)
    * Polygon manipulation
    * ...
* optimization:
  Main functionality for stability analysis and constraint solving
    * Representation of Equation
    * Representation of ValueSet
    * Function which can either consist of an Equation or ValueSet
    * Optimization class which is used to solve constraints defined in Function,
      e.g., find stable regions using a ValueSet => results stored in
      OptimizationResult
    * Parameters represented as parameter set
* parser:
  Classes to parse equations and parameters from strings
* tools:
  General templates used throughout the project, e.g.,
    * An n-dimensional array ndArray
    * A general ThreadPool
    * A general BranchAndBound algorithm
    * Extensions to std::vector in Vec and Mat

The overall stability analysis algorithm is implemented the file glob_stab.cpp
which is compiled to the binary `glob_stab`.

Building glob_stab
==================
Up to now the code has only been compiled and run on Linux. There are only a few
details that would have to be changed to allow compilation on Windows. The main
issue in this context is that it is much easier to install the dependencies,
mainly CGAL, on Linux. Therefore, the software has not been ported to Windows
yet.

Dependencies
------------
* [cmake](https://cmake.org/)
* [g++](https://gcc.gnu.org/)
* [boost](http://www.boost.org/) (tested: 1.58 on Ubuntu 16.04),
  on Ubuntu install
    * libboost-dev
    * libboost-program-options-dev
    * libboost-thread-dev
    * libboost-test-dev
* [CGAL](http://www.cgal.org/) (tested: 4.7-4 on Ubuntu 16.04),
  on Ubuntu install
    * libcgal-dev
* [Doxygen](http://www.doxygen.org) to build the documentatiion

Compilation
-----------
The project uses CMake as build system. For convenience two bash-scripts can be
used to compile the program as debug-build using `build_debug` or as optimized
release-with-debug-symbols build using `build_release`.

In both cases the build is performed in a subdirectory of `build` and the
libraries are copied to `lib` while the executables are copied to `bin` and the
documentation to `doc`.

This behavior may be changed in the top-level `CMakeLists.txt`.
Please note the importance of the compile-time flag `-frounding-math` which is
set in the root-CMakeList.txt and is necessary to allow for a switching of the
rounding mode of the floating-point unit to correctly evaluate interval
arithmetics in Iv.


Using glob_stab
===============
After compiling the program the executable `glob_stab` is used to run a
stability check. The usage of this command is as follows (output of
`glob_stab -h`):

~~~
Help:
  -h [ --help ]         help message

Generic Options:
  -t [ --threads ] arg (=4)             maximum number of parallel threads
  -p [ --parameters ] arg (=cfg/FaFa.p) parameter configuration (file or
                                        inline)
  -c [ --chEq ] arg (=cfg/FaFa.chEq)    characteristic equation (file or
                                        inline)
  -d [ --delta ] arg (=0)               required system damping
  -o [ --output ] arg                   file to store results in
  -a [ --asy ]                          store results as asy (instead of
                                        matlab)
  -s [ --stability ] arg (=1)           stability check (0:none, 1:standard
                                        (boundary), 2:standard (right-half
                                        plane), 3:value-set only (boundary),
                                        4:value-set only (omega+sigma))
  -j [ --complex ]                      complex equation (set s=sigma+j*omega)
  -i [ --no-inverse-valueset ]          do not use inverse to examine unbounded
                                        omega/sigma (with s1/s2 and
                                        finalization)
  -f [ --finalize ]                     for -s1, perform finalizing step
  --omega_max arg (=1000)               Upper bound for omega in value-set only
                                        approach
  --omega_min arg (=-0.0001)            Lower bound for omega in value-set only
                                        approach
  --omega_eps arg (=1e-10)              Epsilon value for omega in
                                        standard-approach
  --sigma_max arg (=1000)               Upper bound for sigma in value-set only
                                        approach
  --sigma_min arg (=0)                  Lower bound for sigma in value-set only
                                        approach
  --sigma_eps arg (=1e-10)              Epsilon value for sigma in
                                        standard-approach
  -v [ --verbosity ] arg (=0)           verbosity of output

Taylor-Model Options:
  --tm_order arg (=10)  order of taylor model
  --tm_center           use domain [-1,1] instead of [0,1]
  --tm_b_off            turn off bernstein bounding
~~~

Robust Stability Analysis
-------------------------
Generally the program is executed by calling

    glob_stab -p <parameter string or file> -c <characteristic equation string or file> -f

which performs an asymptotic stability analysis in two steps (boundary mapping,
stability check for disjoint regions) using a transformation of the imaginary
axis and right half plane which maps this unbounded line/rectangle to two
bounded.

In most cases, additional parameters are necessary:

* `-j`/`--complex`: define \em s as laplace-operator, necessary for time-delay
  systems
* `-o`/`--output <output filename>`: store results in a file instead of printing
  results to stdout in a human-readaby format.
* `-a`/`--asy`: if set then the output is stored as [Asymptote](http://asymptote.sourceforge.net)
  file (only works exactly two parameters with desired resolution set, i.e., 2d
  plots. If not, store results in matlab-format.
* omitt `-f`: For some problems the stability check of disjoint regions takes
  extremely long. By omitting \em -f the program only maps boundaries and
  doesn't perform a stability check of disjoint regions.

Marginal Stability
------------------
Before evaluating the value set and mapping the imaginary axis to the
right half plane the real-root boundary is evaluated (the coefficient a_0) to
check whether the system is marginally stable. Marginal stability can only be
determined using this approach if the coefficient a_0(s=0)==0, i.e., if there is
a root at the origin for s=0. In this case, we first check whether a_1(s=0)==0
as this would imply a double root which means the system is unstable. If
a_1(s=0)!=0 we conclude that there is one root at s=0 and the system is at most
marginally stable.

In this case a small region around zero is excluded for the following boundary
mapping and stability check of disjoint regions. It must be noted that this is
not "safe" as a root might "slip through" this small gap in the imaginary axis
we are not considering (for boundary mapping) or lie in the right-half plane
very close to zero (when checking stability of disjoint regions).

Constraint Solving / Different "Flavours" for Stability Analysis
----------------------------------------------------------------
In addition it is possible to select a slightly different stability analysis
approach or to simply solve a set of inequality constraints by using the
parameter `-s`, e.g.,

* `-s0`: solve general inequality constraints where each equation must be >0.
  Example:

      ./glob_stab -s0 -c 'chEq[0]=y-x*x*x*x; chEq[1]=5-y+x' -p 'x(-4,5,0.01); y(-2,9.1,0.01);' -ao constraint.asy; asy -f png -render 2 constraint.asy

  The result of this command is the following png-image where the resolution in
  x- and y-direction is 1/128th of the interval width due to the specified value
  of 0.01. Within the white region the two specified constraints are positive.

  ![Constraint solving](../images/constraint.png)

  A high-resolution image for the same example where the two resolution-values
  were set to 0.001 resulting in a resolution of 1/1024th of the complete
  interval width can be seen [here](../images/constraint_hq.png).

* `-s1`: default value, first evaluates RRB and checks for marginal stability.
  Then maps the boundary to the parameter space. This is done once for the
  original characteristic equation and once for a transformed characteristic
  equation to completely cover the imaginary axis. If `-f` is specified the
  disjoint regions are then checked for stability using one parameter in each
  disjoint region, again using a transformation to cover the complete right-half
  plane.
* `-s2`: First evaluates RRB and checks for marginal stability. Then maps the
  right-half plane to the parameter space. This is done once for the
  original characteristic equation and once for a transformed characteristic
  equation to completely cover the right-half plane. Usually this takes much
  longer than `-s1` or `-s1 -f` and is not feasible in practice for many
  problems.
* `-s3`: Only evaluates the value set and does not check for marginal stability.
  First maps the boundary to the parameter set and then checks stability for
  each disjoint region if `-f` is also given. In this case the transformation of
  the characteristic function is not performed and we instead check the range
  given by `--omega_min`, `--omega_max` (and `--sigma_min`, `--sigma_max` if
  `-f` is set). Can be useful for cases where convergence is bad for the
  transformed value set.
* `-s4`: Only evaluates the value set and does not check for marginal stability.
  Maps the complete right-half plane to the parameter set. In this case the
  transformation of the characteristic function is not performed and we instead
  check the range given by `--omega_min`, `--omega_max`, `--sigma_min`,
  `--sigma_max`. Can be useful for cases where convergence is bad for the
  transformed value set. Usually you would prefere to use `-s3 -f` as it is
  generally much faster to perform this stability check in two steps.

Specifying the Parameter Set
============================
The parameter set may either be specified by passing a string of parameters or a
filename as argument `-p`. The syntax in both cases is identical:

    a(7.1)
    b1(0.1,1.2);
    AbC_123(-1e-3,2.0,1e-2)

In this case
* **a**: double-parameter '0'
* **b1**: interval [0.1,1.2]
* **AbC_123**: interval [-1e-3,2.0] with a desired resolution of 1e-2

The parameter name must start with an alphabetic character (may be capitalized)
or an underscore and may be followed by an arbitrary number of alphabetic
characters, numbers, or underscores.

A new line or a semicolon or both can be used to separate parameters.

\warning Do not specify **s**, **sigma**, or **omega** as parameters as they are
reserved for internal use!

Specifying characteristic equation
==================================
As for the parameter set the characteristic equation may be passed as inline
string or a filename may be passed as argument `-c`. The systax in both cases is
identical. Here's an example

~~~
c = 3+7*b1
T_d=1e-2
z=exp(-s*T_d)
chEq[0]=3*a+1e-3*b1-AbC_123
chEq[1]=c*c
chEq[2]=2*c*z
~~~

In the end the vector chEq is always used for evaluation. In case of stability
analysis (i.e., the parameter `-s` is set to a positive value) the equation
given above would correspond to the following characteristic equation

\f[
 (3a+1 10^{-3} b1 - AbC_{123} + c^2 \mathbf{s} + 2 c e^{-T_d s} \mathbf{s}^2 = 0
\f]

where \f$c=3+7 b1\f$, i.e., stability for this characteristic equation is
evaluated. In this case the parameter `-j` must be set to make the laplace
operator **s** available which is necessary to specify a time-delay system using
the term \f$z=exp(-s*T_d)\f$.

Note that for `-s0` the program simply solves a number of inequality constraints
and finds the region for which all chEq[i] are positive (obviously this cannot
be evaluated for the equation given above as **s** makes no sense in that case!

Examples
========
Several examples can be found in the directory `examples`. For some of these
examples bash scripts are provided which can be used to run the examples and
generate results.

Copyright and License
=====================
Copyright (C) 2017 Thomas Schauss

    glob_stab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glob_stab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glob_stab. If not, see <http://www.gnu.org/licenses/>.
