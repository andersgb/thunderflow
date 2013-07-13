# coding=UTF-8
## wooot

from collections import namedtuple
import numpy as np
from cpp import pythulip

N = 100 #cells
length = 5.0
dx = length/N
dt = 0.00001

p = pythulip.pythulip('pr_c4', ['nitrogen'])
pythulip.set_SI_units(p)
p.calc()

def temperature_eos_ig(rho, e):
    Cv = 750 #J/(kg K)
    return e/Cv

def pressure_eos_ig(rho, e):
    temperature = temperature_eos_ig(rho, e)
    molarmass = 0.030 # kg/mol
    R = 8.3145 # J / (K mol)
    return rho*R*temperature / molarmass

def eos_update(q):
    rho = q[0]
    vol = dx*1
    velocity = q[1] / rho
    E_internal = (q[2] - 0.5*rho*pow(velocity,2))*vol
    p['var_v'] = vol
    if p.t.number_of_elements('var_n') > 1:
        raise Exception('Not generalized to multicomponent')
    p['var_n'] = rho * p['var_v'][0] / p['cape_molecular_weights'][0]
    p.calc()
    pythulip.iterate_x(p, 'state_u', E_internal)

def pressure_eos(rho, e):
    vol = dx*1
    p['var_v'] = vol
    if p.t.number_of_elements('var_n') > 1:
        raise Exception('Not generalized to multicomponent')
    p['var_n'] = rho * p['var_v'][0] / p['cape_molecular_weights'][0]
    p.calc()
    E = e*rho*vol
    pythulip.iterate_x(p, 'state_u', E)
    return p['state_p'][0]

def calcF(q):
    rho = q[0]
    velocity = q[1] / rho
    e_internal = q[2]/rho-0.5*pow(velocity,2)
    pressure = pressure_eos(rho, e_internal)
    return np.array([rho*velocity, rho*pow(velocity,2)+pressure, (q[2]+pressure)*velocity])

def FORCE(ql, qr):
    fl = calcF(ql)
    fr = calcF(qr)
    q_half = 0.5 * (ql+qr-(dt/dx)*(fr-fl))
    fm = calcF(q_half)
    f_force = 0.25*(fl+2*fm+fr-(dx/dt)*(qr-ql))
    return f_force

def boundary_f(q):
    #velocity zero, but pressure equal
    rho = q[0]
    velocity = q[1] / rho
    e_internal = q[2]/rho-0.5*pow(velocity,2)
    pressure = pressure_eos(rho, e_internal)
    return np.array([0.0, pressure, 0.0])

def musta_evolve(ql, qr):
    fl = calcF(ql)
    fr = calcF(qr)
    f_half = FORCE(ql, qr)
    ql_new = ql - (dt/dx)*(f_half-fl)
    qr_new = qr - (dt/dx)*(fr-f_half)
    return ql_new, qr_new

def evolve(qs, substages=0):
    def find_f(ql_0, qr_0):
        ql = ql_0
        qr = qr_0
        for _ in range(substages):
            ql, qr = musta_evolve(ql, qr)
        return FORCE(ql,qr)
    f_minus = boundary_f(qs[0])
    f_plus = find_f(qs[0], qs[1])
    qs[0] += (dt/dx)*(f_minus - f_plus)
    for i in range(1,N-1):
        f_minus = f_plus
        f_plus = find_f(qs[i], qs[i+1])
        qs[i] += (dt/dx)*(f_minus - f_plus)
    
    f_minus = f_plus
    f_plus = boundary_f(qs[N-1])
    qs[N-1] += (dt/dx)*(f_minus - f_plus)

def initialize():
    p['var_t'] = 300 # K
    p['var_v'] = dx*1 # dx*1m² (arbitrary chosen cross section
    rho = 1.2 #kg/m³
    p['var_n'] = rho * p['var_v'][0] / p['cape_molecular_weights'][0]
    p.calc()
    e_init = p['state_u'][0]
    qs = [np.array([rho, 0.0, e_init*rho]) for _ in range(N)]
    rho_high = 12
    p['var_n'] = rho_high * p['var_v'][0] / p['cape_molecular_weights'][0]
    p.calc()
    qs[-1] = np.array([rho_high, 0.0, p['state_u'][0]*rho_high]) #high pressure in two last CVs
    qs[-2] = np.array([rho_high, 0.0, p['state_u'][0]*rho_high])
    return qs

def sim(step_span):
    qs = initialize()
    while 1:
        yield np.copy(qs)
        for _ in xrange(step_span): evolve(qs)
