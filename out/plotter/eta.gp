set tit "'HTL' vs 'EFF mass'; fix coupling, 1-fnc basis"

fac(Nf) = sqrt(1+Nf/6.)



etaNLL_0(g) = 27.126/log(2.765/g)
etaNLL_2(g) = 86.47/log(2.954/(fac(2)*g))
etaNLL_3(g) = 106.66/log(2.957/(fac(3)*g))
etaNLL_6(g) = 147.63/log(2.994/(fac(6)*g))


set xlab "g"
set ylab "eta *g^4/T^3"

set xr [8e-4:1e2]
set yr [3:1e2]
set log xy
set gr
set key b


set label "n_f = 0" at 4e-2,12
set label "n_f = 3" at 4e-2,20
set label "n_f = 6" at 4e-2,65


p "out/data/etaT3_HTL_nf0_fixed.dat" us 1:($1**4 * $2) w lp lt 1 t "HTL",\
  "out/data/etaT3_HTL_nf3_fixed.dat" us 1:($1**4 * $2) w lp lt 1 t "",\
  "out/data/etaT3_HTL_nf6_fixed.dat" us 1:($1**4 * $2) w lp lt 1 t "",\
  "out/data/etaT3_kappa0.50_nf0_fixed.dat" us 1:($1**4 *$2)   w lp pt 3 lt 3 t "EFF, \153=0.50",\
  "out/data/etaT3_kappa0.50_nf3_fixed.dat" us 1:($1**4 *$2)   w lp pt 3 lt 3 t "",\
  "out/data/etaT3_kappa0.50_nf6_fixed.dat" us 1:($1**4 *$2)   w lp pt 3 lt 3 t "",\
  etaNLL_0(x) lt 4 lw 3 t "NLL",\
  etaNLL_6(x) lt 4 lw 3 t "",\
  etaNLL_3(x) lt 4 lw 3 t ""
pause -1

