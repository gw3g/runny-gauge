etaAMY_Nf0(g) = 27.126/(g**4*log(2.765/g))
etaAMY_Nf2(g) = 86.47/(g**4*log(2.954/(sqrt(4/3)*g)))
etaAMY_Nf3(g) = 106.66/(g**4*log(2.957/(sqrt(3/2)*g)))
etaAMY_Nf6(g) = 147.63/(g**4*log(2.940/(sqrt(2)*g)))

set xl 'alf'
set yl "eta / eta_{NLL}"

set yr [0:2]
set log x

set grid
set key l b

set tit "fix-coupling + eff mass (for u+t denom's) [1-prm Ansatz]"
p "out/data/etaT3_kappa0.25_nf0_fixed.dat"  u 1:($2/etaAMY_Nf0($1)) w lp lt 7 pt 7 t "kappa=0.25, Nf=0",\
  "out/data/etaT3_kappa0.25_nf2_fixed.dat"  u 1:($2/etaAMY_Nf2($1)) w lp lt 6 pt 7 t "kappa=0.25, Nf=2",\
  "out/data/etaT3_kappa0.25_nf6_fixed.dat"  u 1:($2/etaAMY_Nf6($1)) w lp lt 2 pt 7 t "kappa=0.25, Nf=6",\
  "out/data/etaT3_HTL_nf0_fixed.dat"  u 1:($2/etaAMY_Nf0($1)) w lp lt 7 pt 4 t "HTL, Nf=0",\
  "out/data/etaT3_HTL_nf2_fixed.dat"  u 1:($2/etaAMY_Nf2($1)) w lp lt 6 pt 4 t "HTL, Nf=2",\
  "out/data/etaT3_HTL_nf6_fixed.dat"  u 1:($2/etaAMY_Nf6($1)) w lp lt 2 pt 4 t "HTL, Nf=6"
pause -1

