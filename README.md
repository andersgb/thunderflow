A very simple 1D CFD tool that solves the compressible euler equations


It is using a Multi-Stage "MUSTA" solution scheme, see [1] for
details.

The equations are formulated such we simply solve:
```
dq/dt + df(q)/dx = 0
```
where q and f(q) are the vectors:
```
q = [rho, rho*v, E]
f(q) = [rho*v, rho*v^2, (E+p)v]
```

It was developed as a test/learning/demonstration project to see how
small a CFD program in python can be, also demonstrating usage of a
thermodynamic library for equation of state calculations. In this
form here, the thermodynamic library usage is removed and a simple
ideal gas equation of state is used.

TODO:
 - Add multi-component support
 - Extend to two-phase, e.g. like in [2]
 - Find more insanely cool ways to visualize simulations

[1] Titarev, V. A., and E. F. Toro. "MUSTA schemes for
multi-dimensional hyperbolic systems: analysis and improvements."
International journal for numerical methods in fluids 49.2 (2005):
117-148.

[2] Munkejord, Svend Tollak, Steinar Evje, and Tore Flåtten. "A MUSTA
scheme for a nonconservative two-fluid model." SIAM Journal on
Scientific Computing 31.4 (2009): 2587-2622.


Usage instructions:
```
python run.py
python view.py 
```