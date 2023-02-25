import math as m
import sys
import subprocess

def parallel_friction(visc,bondlength,b,L):

    epsilon = 1./m.log(L/(2*b))
    return 2*m.pi*visc*bondlength*epsilon
N = int(sys.argv[1])
L = float(sys.argv[2])
run = int(sys.argv[3])

viscosity = 0.89e-3 # pN us/(nm^2)
rhydro = 3 # nm
kbT = 4.114 # pN nm

bondlength = L/(N-1) # nm

zeta_para = parallel_friction(viscosity,bondlength,rhydro,L) # pN us/nm
zeta_perp = 2*zeta_para # pN us/nm

tau = bondlength**2*zeta_para/kbT #us

dt = 1e-3*tau
TL = zeta_perp/bondlength*L*L*L/kbT/72.0

equilsteps = 3*int(TL/dt)

runsteps = 10*int(TL/dt)
if run <= 10:
    dumpsteps = int(0.01*runsteps)
else:
    dumpsteps = runsteps


print(str(bondlength), str(zeta_para),str(zeta_perp),str(dt),str(equilsteps),
      str(runsteps),str(dumpsteps))

