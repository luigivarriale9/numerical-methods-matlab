# adams_multistep

Adams linear multistep methods in MATLAB: Adams–Bashforth `AdamsB1`–`AdamsB4`
(orders 1–4), Adams–Moulton `AdamsMoulton1`–`AdamsMoulton4` (orders 2–5, implicit
equation solved by fixed-point iteration) and the predictor–corrector `AdamsPECE`
(AB4 predictor, AM4 corrector, order 5). English version of the code I wrote for the
Scientific Computing course, with a few corrections to the original.

`kepler_multistep_convergence.m` measures the observed orders on the Kepler problem
(e = 0.7): position error at t = 3 against the exact solution of Kepler's equation
(`../symplectic_integrators/kepler_position.m`), step counts n = 250 … 4000, the same
as in `kepler_convergence.m`. Starting values come from `ode45` with tight tolerances.

Observed orders: AB1–AB4 0.93, 1.99, 2.88, 3.95; AM1–AM4 2.00, 2.91, 3.97, 5.30;
PECE 5.00. The Python implementation in
[kepler-python](https://github.com/luigivarriale9/kepler-python) gives the same errors:
41 of the 45 values (9 methods × 5 step counts) agree to four significant digits. The
other four are at the finest steps (AM3, AM4 and PECE): three differ in the last digit,
and AM4 at n = 4000 differs by about 3%, where the error (4e-12) is close to the
accuracy of the starting values (DOP853 in Python, `ode45` in MATLAB).
