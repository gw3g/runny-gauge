etaAMY_Nf0(g) = 27.126/(g**4*log(2.765/g))
etaAMY_Nf2(g) = 86.47/(g**4*log(2.954/(sqrt(4/3)*g)))
etaAMY_Nf3(g) = 106.66/(g**4*log(2.957/(sqrt(3/2)*g)))
etaAMY_Nf6(g) = 147.63/(g**4*log(2.940/(sqrt(2)*g)))

set xl 'alf'
set yl "eta / eta_{NLL}"

set yr [0:3]
set log x

set grid
set key l b

set tit "fix-coupling + eff mass (for u+t denom's) [1-prm Ansatz]"
p 'Nf=0, kappa=1.0.dat'  u 1:($2/etaAMY_Nf0(sqrt(4*pi*$1))) w lp lt 7 pt 7 t "kappa=1.0, Nf=0",\
  'Nf=2, kappa=1.0.dat'  u 1:($2/etaAMY_Nf2(sqrt(4*pi*$1))) w lp lt 6 pt 7 t "kappa=1.0, Nf=2",\
  'Nf=6, kappa=1.0.dat'  u 1:($2/etaAMY_Nf6(sqrt(4*pi*$1))) w lp lt 2 pt 7 t "kappa=1.0, Nf=6",\
  'Nf=0, kappa=0.2.dat'  u 1:($2/etaAMY_Nf0(sqrt(4*pi*$1))) w lp lt 7 pt 4 t "kappa=0.2, Nf=0",\
  'Nf=2, kappa=0.2.dat'  u 1:($2/etaAMY_Nf2(sqrt(4*pi*$1))) w lp lt 6 pt 4 t "",\
  'Nf=6, kappa=0.2.dat'  u 1:($2/etaAMY_Nf6(sqrt(4*pi*$1))) w lp lt 2 pt 4 t "",\
  'Nf=0, kappa=0.2, MC=1e6.dat'  u 1:($2/etaAMY_Nf0(sqrt(4*pi*$1))) w lp lt 7 lw 3 pt 4 t "kappa=0.2, Nf=0, MC++",\
  'Nf=2, kappa=1.0_GJ.dat'  u 1:($2/etaAMY_Nf2(sqrt(4*pi*$1))) w lp lt 1 pt 7 t "GJ, Nf=2",\
  'Nf=6, kappa=1.0_GJ.dat'  u 1:($2/etaAMY_Nf6(sqrt(4*pi*$1))) w lp lt 1 pt 7 t "GJ, Nf=6"
pause -1

