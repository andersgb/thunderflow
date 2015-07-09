# coding=UTF-8
## wooot

from collections import namedtuple
import numpy as np

N = 100 #cells
length = 5.0
dx = length/N
dt = 0.00001

#Use ideal gas eos and constant Cv:
# e = T/Cv # e is specific internal energy
# E = rho*e + 0.5*rho*v*v # total energy including kinetic

Cv = 750 #J/(kg K)
molarmass = 0.030 # kg/mol. Air is ca 30 g/mol
R = 8.3145 # J / (K mol)

def temperature_eos_ig(rho, e):
    return e/Cv

def pressure_eos_ig(rho, e):
    temperature = temperature_eos_ig(rho, e)
    return rho*R*temperature / molarmass

def calcF(q):
    """ Calculate the f(q) vector """
    rho = q[0]
    velocity = q[1] / rho
    e_internal = q[2]/rho-0.5*pow(velocity,2)
    pressure = pressure_eos_ig(rho, e_internal)
    return np.array([rho*velocity, rho*pow(velocity,2)+pressure, (q[2]+pressure)*velocity])

def FORCE(ql, qr):
    fl = calcF(ql)
    fr = calcF(qr)
    q_half = 0.5 * (ql+qr-(dt/dx)*(fr-fl))
    fm = calcF(q_half)
    f_force = 0.25*(fl+2*fm+fr-(dx/dt)*(qr-ql))
    return f_force

def boundary_f(q):
    """ Given a solved q inside the domain, return a q with boundary conditions """
    #velocity zero, but pressure equal
    rho = q[0]
    velocity = q[1] / rho
    e_internal = q[2]/rho-0.5*pow(velocity,2)
    pressure = pressure_eos_ig(rho, e_internal)
    return np.array([0.0, pressure, 0.0])

def musta_evolve(ql, qr):
    fl = calcF(ql)
    fr = calcF(qr)
    f_half = FORCE(ql, qr)
    ql_new = ql - (dt/dx)*(f_half-fl)
    qr_new = qr - (dt/dx)*(fr-f_half)
    return ql_new, qr_new

def evolve(qs, substages=0):
    """ Do a single timestep """
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

def init_example1():
    """Set an initial state where the last two CVs have a high pressure"""
    T = 300 # K
    rho = 1.2 #kg/m³
    E = rho*Cv*T
    qs = [np.array([rho, 0.0, E]) for _ in range(N)]

    #high pressure in two last CVs:
    rho_high = 12
    E_high = rho_high*Cv*T
    qs[-1] = np.array([rho_high, 0.0, E_high])
    qs[-2] = np.array([rho_high, 0.0, E_high])
    return qs

def sim(step_span):
    #each control volume solves for q = [rho, rho*v, E]

    #qs is a vector containing over all control volumes, with one q
    #for each. I.e. qs[CV_INDEX][VAR_INDEX]
    qs = init_example1()
    while 1:
        yield np.copy(qs)
        for _ in xrange(step_span): evolve(qs)
