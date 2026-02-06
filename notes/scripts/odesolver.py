import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif"
})

# Physical constants (you may change them)
K = 1.0
G = 1.0
rho_c = 1.0
n = 1
A = K * (1 + 1/n) * rho_c**(1/n - 1) / (4*np.pi*G)

def lane_emden_general(r, y, n):
    """
    Convert equation to first order system.
    y = [rho, drho/dr]
    """
    phi, dphi = y

    if phi <= 0:
        return [0, 0]   # stop integration when density becomes non-positive

    # coefficient A(phi)

    # Compute second derivative from original equation:
    # (1/r^2) d/dr (r^2 A dphi/dr) = -phi^n
    term = -(phi**n) - (2/r) * dphi
    d2phi = term

    return [dphi, d2phi]


def solve_for_n(n, r_max=10, phi0=1.0):
    """
    Solve the equation for a given n.
    """
    # Initial conditions:
    # rho(0) = rho0
    # drho/dr (0) = 0 but numerically we start at a small r
    r0 = 1e-6
    phi_init = phi0
    dphi_init = 0.0

    sol = solve_ivp(lambda r, y: lane_emden_general(r, y, n),
                    [r0, r_max],
                    [phi_init, dphi_init],
                    max_step=0.01,
                    dense_output=True,
                    events=lambda r, y: y[0])  # stop at rho = 0

    return sol


# --- Run for several n values ---
n_values = [1e-6, 1, 2, 3, 4, 5]

plt.figure(figsize=(8, 6))

for n in n_values:
    sol = solve_for_n(n)
    r = sol.t
    phi = sol.y[0]

    plt.plot(r, phi, label=f"n = {n}")

plt.ylim(0, 1.05)
plt.xlabel(r"$\eta$", fontsize = 14)
plt.ylabel(r"$\phi(\eta)$", fontsize = 14)
plt.title("Solutions of generalized Lane-Emden equation", fontsize = 14)
plt.legend()
plt.grid(True)
plt.savefig('lane_enden.pdf', transparent = True)
plt.show()
