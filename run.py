import numpy as np

import BBGKY as method
import model.spinboson_1D as model


TIERS = 20
FILTER_THRESHOLD = 1e-12
FILTER_INTERVAL = 50


par = model.parameters
ndof = par.ndof
NSteps = par.NSteps
NStates = par.NStates
dt = par.dt
nskip = par.nskip
initState = getattr(par, "initState", 0)


N, up1, up2, dn1, dn2, m1_arr, m2_arr = method.build_ado_structure(
    ndof, tier=TIERS
)
print(f"nmax = {N}")


omega = par.omega

par.Hs = model.Hs()
par.Qs = model.Qs()
par.up1, par.up2 = up1, up2
par.dn1, par.dn2 = dn1, dn2
par.m1_arr, par.m2_arr = m1_arr, m2_arr
par.gamma = np.sum(
    (m1_arr + m2_arr) * np.imag(omega)
    - 1j * np.real(omega) * (m1_arr - m2_arr),
    axis=1,
)


rhot = np.zeros((NStates, NStates, N), dtype=complex)
rhot[initState, initState, 0] = 1.0


active = None
n_active = N

with open("Pt.txt", "w") as f:
    for k in range(NSteps):
        if k % nskip == 0:
            pct = 100.0 * k / NSteps
            print(f"\r  {pct:5.1f}%  t={k * dt:.3f}  nADOs={n_active}/{N}", end="", flush=True)
            f.write(f"{round(k * dt, 3)}\t")
            for i in range(NStates):
                f.write(f"{np.real(rhot[i, i, 0])}\t")
            f.write("\n")

        if FILTER_INTERVAL > 0 and k % FILTER_INTERVAL == 0:
            active = method.filter_ados(rhot, FILTER_THRESHOLD)
            n_active = int(active.sum())

        rhot = method.prop(rhot, par, active=active)

print(f"\r  100.0%  t={NSteps * dt:.3f}  nADOs={n_active}/{N}  - done")
