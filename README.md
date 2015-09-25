# (shear) viscosity in quark-gluon plasma

transport coefficients in hot gauge theories, eta (T,g) with Nf = 0, ... 6

_Usage_: (produces CSV data files, default 20 points)
* ```eval_g(gmin, gmax)``` evaluates \eta(g)/T^3
* ```eval_T(Tmin, Tmax)``` evaluates \eta(T)/T^3


leading-log [AMY1](http://arxiv.org/abs/hep-ph/0010177)
& beyond [AMY2](http://arxiv.org/abs/hep-ph/0302165)


  file          |   purpose
----------------|-------------:
  main.c        |   unexpanded eta calc             [MAIN FILE]
  core.h        |   header file
  thermal.c     |   quantum statistics
  basis.c       |   expand chi(p) in basis
  amy_prepInt.c |   organise integration (q,omega)
  stu_prepInt.c |                        (s, t)
  operators.c   |   define functional operators
  qcd.c         |   matrix elements
  htl.c         |   hard thermal loops

A *rough* UML diagram is shown in "~/doc/layout.jpg".

**monte carlo** integrator using gsl implementation of [VISER](https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS)

OUTPUT: to "~/out/", use GLE to plot.

