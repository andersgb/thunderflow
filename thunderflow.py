## wooot

from collections import namedtuple
import numpy as np

N = 10 #cells
length = 5.0
dx = length/N
dt = 0.00001

Q = namedtuple('Q', 'rho rho_u')
F = namedtuple('F', 'rho_u rho_usq_p')

def pressure_eos(rho):
    temperature = 300
    molarmass = 0.030 # kg/mol
    R = 8.3145 # J / (K mol)
    return rho*R*temperature / molarmass

def calcF(q):
    rho = q[0]
    velocity = q[1] / rho
    pressure = pressure_eos(rho)
    return np.array([rho*velocity, rho*velocity*velocity+pressure])

def FORCE(ql, qr):
    fl = calcF(ql)
    fr = calcF(qr)
    q_half = 0.5 * (ql+qr-(dt/dx)*(fr-fl))
    fm = calcF(q_half)
    f_force = 0.25*(fl+2*fm+fr-(dx/dt)*(qr-ql))
    return f_force

def boundary_f(q):
    #velocity zero, but pressure equal
    pressure = pressure_eos(q[0])
    return np.array([0.0, pressure])

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
    qs = [np.array([1.2, 0.0]) for _ in range(N)]
    qs[-1] = np.array([12, 0.0])
    qs[-2] = np.array([12, 0.0])
    return qs


