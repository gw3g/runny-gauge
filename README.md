# Shear viscosity in the quark-gluon plasma

The shear viscosity for QCD, η(T),
is revisited with particular focus on the relevant scale(s) that determine
α, the strong-coupling.
Perturbation theory was previously deemed incompatible with the inference from
heavy ion collisions (which indicated "perfect fluid" characteristics) for
reasons we explain in
[1704.06284](https://arxiv.org/abs/1704.06284).
Also in the aforementioned preprint, we give the revised temperature dependence
for T≈T<sub>c</sub>, which is compatible with other estimates from lattice QCD
and hydrodynamic simulations.

** Background / kinetic theory **

The QCD Boltzmann equation (in the quenched case) was first linearised by
[Baym et al.](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.64.1867)
who used hard thermal loop (HTL) insertions to screen the collision term.
Their calculation fixed the overall *leading-log* prefactor, and has been
extended to include n<sub>f</sub>>0 flavours of massless quarks
[plus other gauge theories, besides SU(N)] due to Arnold, Moore & Jaffe (AMY)
  [[Part I]](http://arxiv.org/abs/hep-ph/0010177).
AMY also showed how to go beyond logarithmic accuracy 
 by incorporating inelastic processes in
 [[Part II]](http://arxiv.org/abs/hep-ph/0302165).
 This settled the LO result for η, as an expansion in α and allowed the
 *next-to-leading-log* coefficients to be calculated.


## Data

**fixed-α calculation**

Here we present results for the LO viscosity, which are accurate to 1% (we only dropped inelastic contributions
and used a single-function _ansatz_ for the off equilibrium solution).
Data is saved under ``out/data/``, where file names
indicate etaT3 (for η/t<sup>3</sup>).
Firstly, the files tagged 'fixed' give the parametric dependence on the coupling
parameter g=(4πα)<sup>1/2</sup>. Some useful links are given below:

* [etaT3_HTL_nf0_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf0_fixed.dat)
* [etaT3_kappa0.50_nf0_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf0_fixed.dat)
* [etaT3_kappa0.25_nf0_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.25_nf0_fixed.dat)
* [etaT3_HTL_nf2_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf2_fixed.dat)
* [etaT3_kappa0.50_nf2_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf2_fixed.dat)
* [etaT3_kappa0.25_nf2_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.25_nf2_fixed.dat)
* [etaT3_HTL_nf3_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf3_fixed.dat)
* [etaT3_kappa0.50_nf3_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf3_fixed.dat)
* [etaT3_kappa0.25_nf3_fixed.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.25_nf3_fixed.dat)

**With Running**

Here we implement the scheme devised by AMY (minus the inelastic processes),
with a coupling that depends on the virtuality Q<sup>2</sup>=ω<sup>2</sup>-q<sup>2</sup>.
An effective version of the one-loop running is used, to account for both
s-channel (Q<sup>2</sup>>0) and t- or u-channel (Q<sup>2</sup><0) processes.
To account for the singularity at the Landau pole,
a maximal value α≤[1...10] which turns out to be of little importance.

![eta with running](https://github.com/gw3g/runny-gauge/blob/master/out/eta_running.png)

The files tagged as 'running' give the temperature dependence dependence
(T in units of Λ<sub>QCD</sub>).

* [etaT3_HTL_nf0_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf0_running.dat)
* [etaT3_kappa0.50_nf0_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf0_running.dat)
* [etaT3_HTL_nf2_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf2_running.dat)
* [etaT3_kappa0.50_nf2_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf2_running.dat)
* [etaT3_HTL_nf3_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_HTL_nf3_running.dat)
* [etaT3_kappa0.50_nf3_running.dat](https://github.com/gw3g/runny-gauge/blob/master/out/data/etaT3_kappa0.50_nf3_running.dat)

## Usage

I have designed this piece of code for my personal use, simply to calculate η (for QCD).
However, this has become a central quantity in heavy-ion physics and these results may be
of wider interest.

**Visualisation**

Use [gnuplot](http://www.gnuplot.info) to plot data (scripts are under ``out/plotter``).
For example, to show the temperature dependent viscosity for 3-flavour QCD:
```
make 'NF=3' temperature
```

**Details on 5-dimensional integration**

Use the makefile to compile the binary, execute with
```
./bin/eta NF
```
where ```NF```={0,1,...} is the number of light quark flavours.
(The default value is ```Nf=0```.)

* ```eval_g(gmin, gmax)``` evaluates η(g)/T<sup>3</sup>
* ```eval_T(Tmin, Tmax)``` evaluates η(T/Λ)/T<sup>3</sup>
* ```rate_E(Emin, Emax, Temp)``` evaluates absorbtion rate R(E)/T

Using HTL function (set ```HTL=1```) OR eff. mass ```mu^2 = kappa*mD^2``` (set ```HTL=0```)


Integrator using gsl implementation of (pick one).
The absolute value of ``(int) calls`` indicates how many times the integrand
is sampled (the sign determines which codes below are used).

  * [VEGAS](https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html#VEGAS),  ``calls > 0``
  * [Cubature](http://ab-initio.mit.edu/wiki/index.php/Cubature),   ``calls < 0`` ~ *recommended*

Note that ``±1e5`` is _just_ stable for n<sub>f</sub>=0, and should
be higher by a factor ~10<sup>3</sup> with quarks.
