# Numerical methods in MATLAB

MATLAB implementations from my university work: integrators for Hamiltonian
systems (from my M.Sc. thesis in Mathematics). Adams multistep methods,
optimisation and least squares (Scientific Computing course, joint assignments)
will be added. The Python versions with automated tests are in
[kepler-python](https://github.com/luigivarriale9/kepler-python).

## symplectic_integrators

Four one-step methods applied to the Kepler problem (eccentricity 0.7):

| File | Method | Order |
|---|---|---|
| `Euler_Explicit.m` | explicit Euler | 1 |
| `Symplectic_euler.m` | symplectic Euler | 1 |
| `Midpoint_Implicit.m` | implicit midpoint (`fsolve`) | 2 |
| `Stormer_Verlet.m` | Stormer-Verlet | 2 |

`kepler_convergence.m` computes the global error at t = 3 against the exact solution
(`kepler_position.m`, Kepler's equation solved with Newton) and prints the observed order:

```
Explicit Euler     0.93 (expected 1)
Symplectic Euler   0.95 (expected 1)
Implicit midpoint  2.00 (expected 2)
Stormer-Verlet     2.00 (expected 2)
```

   ![Order of convergence](figures/kepler_convergence.png)

Requires the Optimization Toolbox for `fsolve`. The script ends with an `assert`
on each order, so it also works as a check. Tested in GNU Octave 8.4.
