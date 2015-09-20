# (shear) Viscosity in quark-gluon plasma

transport coefficients in hot gauge theories, eta (T,g) with Nf = 0, ... 6

leading-log [AMY1](http://arxiv.org/abs/hep-ph/0010177)
& beyond [AMY2](http://arxiv.org/abs/hep-ph/0302165)

```
  file          :   purpose
-------------------------------------------------------------
  main.c        :   unexpanded eta calc             [MAIN FILE]
  core.h        :   header file
  thermal.c     :   quantum statistics
  basis.c       :   expand chi(p) in basis
~/integrand/
  amy_prepInt.c :   organise integration (q,omega)
  stu_prepInt.c :                        (s, t)
  operators.c   :   define functional operators
~/x-sec/
  qcd.c         :   matrix elements
  htl.c         :   hard thermal loops
```

Currently, the function main() will calculate eta(g) in units
of T for Nf = 0,...,6 . A *rough* UML diagram is shown in 
"~/doc/layout.png".

OUTPUT: to "~/out/", csv files. Use "~/out/plot.nb" to view.

