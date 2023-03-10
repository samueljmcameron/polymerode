import math as m
import sys
import subprocess

def parallel_friction(visc,bondlength,b,L):

    epsilon = 1./m.log(L/(2*b))
    return 2*m.pi*visc*bondlength*epsilon



def set_params(N,bondlength=4.0,viscosity=0.89e-3,rhydro=2.0,kbT=4.114):
    """
    Set timescale-related parameters given input data

    Params
    ------

    N : int
        number of beads in the polymer
    bondlength : double (optional)
        length between beads. default value is 4 (nm), reasonable size for polymer
        with persistence length of 50 (nm)
    viscosity : double (optional)
        solvent viscosity. default value is 8.9e-4 (pN us/nm^2), water at 25 deg C
    rhydro : double (optional)
        hydrodynamic radius of polymer (roughly). default value is 2 (nm), reasonable
        for ds-DNA
    kbT : double (optional)
        temperature. default value is 4.114 (pN nm)
    

    Returns
    -------

    output : dict
        Return parameters as parallel and perpendicular friction,
        short time-scale (usually need a timestep 1e-3 times this), and long
        time-scale (rotational diffusion of rod of same length)
    """
    


    outputs = {}



    L = bondlength*(N-1) # nm

    zeta_para = parallel_friction(viscosity,bondlength,rhydro,L) # pN us/nm
    outputs['friction_para'] = zeta_para
    outputs['friction_perp'] = 2*zeta_para # pN us/nm

    outputs['short_tscale'] = bondlength**2*zeta_para/kbT #us

    outputs['long_tscale'] = 2*zeta_para/bondlength*L*L*L/kbT/72.0

    return outputs
