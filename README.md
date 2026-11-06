# Numerical methods in MATLAB

MATLAB implementations from my university work: integrators for Hamiltonian
systems (from my M.Sc. thesis in Mathematics) and Adams multistep methods
(Scientific Computing course). Optimisation and least squares will be added. The
Python versions with automated tests are in
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
on each order, so it also works as a check.

## adams_multistep

Adams–Bashforth `AdamsB1`–`AdamsB4` (orders 1–4), Adams–Moulton
`AdamsMoulton1`–`AdamsMoulton4` (orders 2–5, implicit equation solved by fixed-point
iteration) and the predictor–corrector `AdamsPECE` (AB4 predictor, AM4 corrector,
order 5), from the Scientific Computing course, with a few corrections.

`kepler_multistep_convergence.m` measures the observed orders on the Kepler problem
(position error at t = 3 against the exact solution, n = 250 … 4000, the same as
for the symplectic integrators):

| Method | AB1 | AB2 | AB3 | AB4 | AM1 | AM2 | AM3 | AM4 | PECE |
|---|---|---|---|---|---|---|---|---|---|
| Observed order | 0.93 | 1.99 | 2.88 | 3.95 | 2.00 | 2.91 | 3.97 | 5.30 | 5.00 |

![Adams methods on the Kepler problem](figures/kepler_multistep_convergence.png)

The errors agree with the Python implementation in
[kepler-python](https://github.com/luigivarriale9/kepler-python): 41 of the 45 values
(9 methods × 5 step counts) match to four significant digits; the other four are at the
finest steps, and AM4 at n = 4000 differs by about 3%, where the error (4e-12) is close
to the accuracy of the starting values (`ode45` here, DOP853 in Python).

Tested in MATLAB R2024b and GNU Octave 8.4.