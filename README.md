# Shear viscosity in the quark-gluon plasma (v1.0)

The transport coefficient η (shear viscosity) for
QCD is revisited with particular focus on the relevant scale(s) that determine
α, the strong-coupling.
Perturbation theory was previously deemed incompatible with the inference from
heavy ion collisions (which indicated "perfect fluid" characteristics) for
reasons we explain in
[1704.06284](https://arxiv.org/abs/1704.06284).
Also in the aforementioned preprint, we give the revised temperature dependence
for T≈T<sub>c</sub>, which is compatible with other estimates from lattice QCD
and hydrodynamic simulations.

## Background / kinetic theory

**fixed-α calculation**

The QCD Boltzmann equation (in the quenched case) was first linearised by
[Baym et al.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.64.1867)
who used hard thermal loop (HTL) insertions to screen the collision term.
Their calculation fixed the overall *leading-log* prefactor, has been
extended to include n<sub>f</sub>>0 flavours of massless quarks
(and for other gauge groups) due to Arnold, Moore & Jaffe (AMY)
  [[Part I]](http://arxiv.org/abs/hep-ph/0010177).
AMY also showed how to go beyond logarithmic accuracy in
 [[Part II]](http://arxiv.org/abs/hep-ph/0302165),
 and incoporated inelastic processes.
 That established the rigourous 'weak coupling' NLL formula.

<!--![eta-over-NLL](https://github.com/gw3g/runny-gauge/blob/master/out/eta-to-NLL.png)-->

**With Running**

Here we implement the scheme devised by AMY (minus the inelastic processes),
with a coupling that depends on the virtuality Q<sup>2</sup>=ω<sup>2</sup>-q<sup>2</sup>.
An effective version of the one-loop running is used, to account for both
s-channel (Q<sup>2</sup>>0) and t- or u-channel (Q<sup>2</sup><0) processes.
To account for the singularity at the Landau pole,
a maximal value α≤[1...10] which turns out to be of little importance.

<!--![eta with running](https://github.com/gw3g/runny-gauge/blob/master/out/eta_running.png)-->

## Usage

Data is saved under ``out/data/``, where file names
indicate etaT3 (for η/T<sup>3</sup>).
Firstly, the files tagged 'fixed' give the parametric dependence on the coupling
parameter g=(4πα)<sup>1/2</sup>.
The files tagged as 'running' give the temperature dependence dependence
(T in units of Λ<sub>QCD</sub>).

* ```eval_g(gmin, gmax)``` evaluates η(g)/T<sup>3</sup>
* ```eval_T(Tmin, Tmax)``` evaluates η(T/Λ)/T<sup>3</sup>

Using HTL function (set ```HTL=1```) OR eff. mass ```mu^2 = kappa*mD^2``` (set ```HTL=0```)

**Visualisation**

Plots are made using [GLE](http://glx.sourceforge.net/)
(scripts are under ``out/plotter``),
see the Makefile.

**Details on 5-dimensional integration**

Integrator using gsl implementation of (pick one).
The absolute value of ``(int) calls`` indicates how many times the integrand
is sampled (the sign determines which codes below are used).

  * [VEGAS](https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS),  ``calls > 0``
  * [Cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature),   ``calls < 0`` ~ *recommended*

Note that ``±1e5`` is _just_ stable for n<sub>f</sub>=0, and should
be higher by a factor ~10<sup>3</sup> with quarks.
