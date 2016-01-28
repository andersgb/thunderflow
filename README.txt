A very simple 1D CFD tool that solves the compressible euler equations

It is using a Multi-Stage "MUSTA" solution scheme, see [1] for
details.

The equations are formulated such we simply solve:
dq/dt + df(q)/dx = 0

where q and f(q) are the vectors:

q = [rho, rho*v, E]
f(q) = [rho*v, rho*v^2, (E+p)v]

TODO:
 - Add multi-component support
 - Extend to two-phase, e.g. like in [2]
 - Find more insanely cool ways to visualize simulations

[1] MUSTA schemes for multi-dimensional hyperbolic systems: analysis
    and improvements. V.A. Titarev and E.F. Toro
[2] A MUSTA scheme for a nonconservative two-fluid model. S.T.
    Munkejord, S. Evje, T. Flåtten

Usage instructions:
python run.py
python view.py 
