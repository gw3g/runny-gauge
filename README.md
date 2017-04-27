# Shear viscosity in the quark-gluon plasma (v1.0)

The transport coefficient η (shear viscosity) for
QCD is revisited with particular focus on the relevant scale that determines
α, the strong-coupling.
Perturbation theory was previously deemed incompatible with the inference from
heavy ion collisions (which indicated "perfect fluid" characteristics) for
reasons we explain in
[1704.06284](https://arxiv.org/abs/1704.06284).
Also in the aforementioned preprint, we give the revised temperature dependence
for T≈T<sub>c</sub>, which is compatible with other estimates from lattice QCD
and hydrodynamic simulations.

----
## usage

_Usage_: (produces DAT data files)
* ```eval_g(gmin, gmax)``` evaluates \eta(g)/T^3
* ```eval_T(Tmin, Tmax)``` evaluates \eta(T)/T^3

Using HTL function (set ```HTL=1```) OR eff. mass ```mu^2 = kappa*mD^2``` (set ```HTL=0```)

Data is saved under ``~/out/data/``
The file names indicate etaT3 (for η/T<sup>3</sup>)

----
## parametric α-dependence

The QCD Boltzmann equation (in the quenched case) was first linearised by
[Baym et al.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.64.1867)
who used hard thermal loop (HTL) insertions to screen the collision term.
Their calculation fixed the overall *leading-log* prefactor, has been
extended to include n<sub>f</sub>>0 flavours of massless quarks
(and for other gauge groups) due to Arnold, Moore & Jaffe (AMY)
  [[Part I]](http://arxiv.org/abs/hep-ph/0010177).
AMY also showed how to go beyond logarithmic accuracy in
 [[Part II]](http://arxiv.org/abs/hep-ph/0302165),
 also accomodating inelastic processes.

 given in terms of g = (4pi )

----
## resummed pert



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
