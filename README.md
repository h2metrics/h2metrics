h2metrics
=========

What is it?
-----------

This code provides tools for Riemannian shape analysis with second order Sobolev metrics on closed plane curves. It can be used to solve the initial and boundary value problems for geodesics and to compute Karcher means. It is able to factor out reparameterizations, translations, rotations and scalings of the curves. It works for both open and closed planar curves.

For details we refer to the our papers

    @article{BBHM2017,
      author  = {Martin Bauer, Martins Bruveris, Philipp Harms, Jakob M{\o}ller-Andersen},
      title   = {A Numerical Framework for {S}obolev Metrics on the Space of Curves},
      journal = {SIAM J. Imaging Sci.},
      year    = {2017},
      volume  = {10},
      number  = {1},
      pages   = {47--73}
    }

    @misc{BBCH2018,
      author  = {Martin Bauer, Martins Bruveris, Nicolas Charon, Jakob M{\o}ller-Andersen},
      title   = {A relaxed approach for curve matching with elastic metrics},
      year    = {2018},
      note    = {Preprint available at arXiv:1803.10893}
    }

Please cite our papers in your work.

Dependencies
------------

* MATLAB Curve Fitting Toolbox
* MATLAB Optimization Toolbox

The code incorporates the following libraries

* Manopt library (www.manopt.org)
* Hanso (https://cs.nyu.edu/overton/software/hanso/)
* export_fig (https://github.com/altmany/export_fig)
* fshapesTk (https://github.com/fshapes/fshapesTk)

The code was tested on MATLAB R2016b.

Usage
-----

See the file "example.m" for an example of how to use the code.

Licence
-------

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program. If not, see http://www.gnu.org/licenses/.

Contacts
--------

* Martin Bauer (martin dot bauer at tuwien dot ac dot at)
* Martins Bruveris (martins dot bruveris at brunel dot ac dot uk)
* Nicolas Charon (charon at cis dot jhu dot edu)
* Philipp Harms (philipp dot harms at stochastik dot uni-freiburg dot de)
* Jakob Møller-Andersen (jakmo at dtu dot dk)
