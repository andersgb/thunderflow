## wooot

from collections import namedtuple
import numpy as np

N = 10 #cells
length = 5.0
dx = length/N
dt = 0.00001

def temperature_eos(rho, e):
    Cv = 750 #J/(kg K)
    return e/Cv

def pressure_eos(rho, e):
    temperature = temperature_eos(rho, e)
    molarmass = 0.030 # kg/mol
    R = 8.3145 # J / (K mol)
    return rho*R*temperature / molarmass

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
    e_init = 750*300
    qs = [np.array([1.2, 0.0, e_init*1.2]) for _ in range(N)]
    qs[-1] = np.array([12, 0.0, e_init*12])
    qs[-2] = np.array([12, 0.0, e_init*12])
    return qs

