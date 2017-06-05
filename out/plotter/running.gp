set tit "'renorm. LO' vs NLL; running alpha, 1-fnc basis"

Nf = ARG1

alpha(Q,Lambda) = 4*pi/((11-2*Nf/3)*2*log(Q/Lambda))

fac(Nf) = sqrt(1+Nf/6.)
s0 = (16+10.5*Nf)*4*pi**2/90

# See (4.17) in hep-ph/0302165
if (Nf==0) {etaNLL(g) = 27.126/(g**4 *log(2.765/g));}
if (Nf==1) {etaNLL(g) = 60.808/(g**4 *log(2.7/(fac(1)*g)));}
if (Nf==2) {etaNLL(g) = 86.47/(g**4 *log(2.954/(fac(2)*g)));}
if (Nf==3) {etaNLL(g) = 106.66/(g**4 *log(2.957/(fac(3)*g)));}
if (Nf==4) {etaNLL(g) = 122.96/(g**4 *log(2.954/(fac(4)*g)));}
if (Nf==5) {etaNLL(g) = 136.38/(g**4 *log(2.947/(fac(5)*g)));}
if (Nf==6) {etaNLL(g) = 147.63/(g**4 *log(2.994/(fac(6)*g)));}

set xl "T/Lambda"
set yl "eta / T^3"

set xr [1:9]
set yr [0:15+15*Nf]
#set log y

set key t l
set gr


set label "n_f = 0" at 4e-2,12

# set tit "vary alf (T=0.2 GeV)"
p etaNLL(sqrt(4*pi*alpha(2*pi*x,0.5))) lt 1 t "Lambda = [0.5,1,2]Tc",\
  etaNLL(sqrt(4*pi*alpha(2*pi*x,1)))   lt 1 t "",\
  etaNLL(sqrt(4*pi*alpha(2*pi*x,2)))   lt 1 t "",\
  "out/data/etaT3_HTL_nf".Nf."_running.dat" u 1:2:3:4 w err lt 7 t " t^* = -[.5,2]T^2",\
  "out/data/etaT3_HTL_nf".Nf."_running.dat" u 1:5 lt 6 t " omni-HTL",\
  s0 lt 5 t "s_0/T^4"
pause -1

