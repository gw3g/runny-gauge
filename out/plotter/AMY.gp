Nf = 0

# Lambda in units of Tc, expect Lambda ~ Tc, cf http://arxiv.org/pdf/hep-ph/0601119.pdf
alpha(Q,Lambda) = 4*pi/((11-2*Nf/3)*2*log(Q/Lambda))

etaAMY(g) = 27.126/(g**4*log(2.765/g))

s0 = 16*4*pi**2/90


set xl 'T'
set yl "eta / s"

set xr [1:4]
set yr [0.04:4]
#set log y

set key b
set grid


# set tit "vary alf (T=0.2 GeV)"
p etaAMY(sqrt(4*pi*alpha(2*pi*x,0.5)))/s0 t "Lambda = 0.5Tc",\
  etaAMY(sqrt(4*pi*alpha(2*pi*x,1)))/s0 t "Lambda = Tc",\
  etaAMY(sqrt(4*pi*alpha(2*pi*x,2)))/s0 t "Lambda = 2Tc",\
  "Boyd1996_pes-spline.dat" u 1:(etaAMY(sqrt(4*pi*alpha(2*pi*$1,0.5)))/$4) t " 0.5",\
  "Boyd1996_pes-spline.dat" u 1:(etaAMY(sqrt(4*pi*alpha(2*pi*$1,1)))/$4) t " 1.0",\
  "Boyd1996_pes-spline.dat" u 1:(etaAMY(sqrt(4*pi*alpha(2*pi*$1,2)))/$4) t " 2.0"

pause -1

