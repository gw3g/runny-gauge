# (shear) viscosity in quark-gluon plasma

transport coefficients in hot gauge theories, eta (T,g) with Nf = 0, ... 6

_Usage_: (produces DAT data files)
* ```eval_g(gmin, gmax)``` evaluates \eta(g)/T^3
* ```eval_T(Tmin, Tmax)``` evaluates \eta(T)/T^3

Using HTL function (set ```HTL=1```) OR eff. mass ```mu^2 = kappa*mD^2``` (set ```HTL=0```)

## literature

leading-log [AMY1](http://arxiv.org/abs/hep-ph/0010177)
& beyond [AMY2](http://arxiv.org/abs/hep-ph/0302165)

## organisation

  file          |   purpose
----------------|-------------:
  main.c        |   unexpanded eta calc             [MAIN FILE]
  core.h        |   header file
  thermal.c     |   quantum statistics
  basis.c       |   expand chi(p) in basis
  integrand.c   |   organise integration (q,omega) / (t,x)
  operators.c   |   define functional operators [^1]
  qcd.c         |   matrix elements
  htl.c         |   hard thermal loops

[^1] **monte carlo** integrator using gsl implementation of (pick one)
* [VEGAS](https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS)
* [Cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature)

## visualisation

data (comma separated) is written to ``~/out/data/``. 
Pretty pictures  using [GLE](http://glx.sourceforge.net/).

