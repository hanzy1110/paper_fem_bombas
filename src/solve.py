import numpy as np
from scipy.special import expit


def piston_action(alfa, p_linea, p_servicio):
    alfa %= 2 * np.pi
    return (p_servicio - p_linea) * expit((alfa - np.pi) * 20)


def build_indep(
    alfa: np.float32,
    beta: np.float32,
    omega: np.float32,
    mb: np.float32,
    mp: np.float32,
    R: np.float32,
    Ib: np.float32,
    Lg: np.float32,
    L: np.float32,
    D_plunger,
    p_linea,
    p_servicio,
) -> np.ndarray:
    p_alfa = piston_action(alfa, p_linea, p_servicio)

    dbdt = (R / L) * (np.cos(alfa) / np.cos(beta)) * omega
    d2bdt2 = (dbdt**2 - omega**2) * np.tan(beta)

    b0 = mb * Lg * (R / L) * omega**2 * np.sin(alfa)
    b1 = (
        -mb
        * R
        * omega**2
        * (
            np.cos(alfa)
            + (1 - Lg / L)
            * (
                (R / L) * (np.cos(alfa) ** 2 / np.cos(beta) ** 3)
                - np.sin(alfa) * np.tan(beta)
            )
        )
    )
    b2 = Ib * d2bdt2
    b3 = 0
    b4 = mp * R * omega**2 * (
        (np.cos(alfa + beta) / np.cos(beta))
        + (R / L) * (np.cos(alfa) ** 2 / np.cos(beta) ** 3)
    ) - p_alfa * (D_plunger**2 * np.pi / 4)

    return np.array([b0, b1, b2, b3, b4])


def build_matrix(L: np.float32, Lg: np.float32, beta: np.float32) -> np.ndarray:
    m = np.eye(5)
    m[4, 4] = 0
    m[4, 3] = 1

    m[3, 3] = 0
    m[3, 4] = -1
    m[3, 2] = 1

    m[0, 2] = 1
    m[1, 3] = 1

    vec = np.array(
        [
            -(L - Lg) * np.sin(beta),
            (L - Lg) * np.sin(beta),
            Lg * np.cos(beta),
            -Lg * np.sin(beta),
            0,
        ]
    )

    m[2] = vec
    return m


def solve_system(
    alfa,
    beta,
    omega,
    mb: np.float32,
    mp: np.float32,
    R: np.float32,
    Ib: np.float32,
    L,
    Lg,
    D_plunger,
    p_linea,
    p_servicio,
):
    A = build_matrix(L, Lg, beta)
    b = build_indep(
        alfa, beta, omega, mb, mp, R, Ib, Lg, L, D_plunger, p_linea, p_servicio
    )

    sol = np.linalg.solve(A, b) / 1000

    R = np.sqrt((sol[0] + sol[2]) ** 2 + (sol[1] + sol[3]) ** 2)
    theta_R = np.arctan((sol[1] + sol[3]) / (sol[0] + sol[2]))

    ret_dict = {
        "module_A": np.sqrt(sol[0] ** 2 + sol[1] ** 2),
        "phase_A": np.arctan2(sol[1], sol[0]),
        "module_B": np.sqrt(sol[2] ** 2 + sol[3] ** 2),
        "phase_B": np.arctan2(sol[3], sol[2]),
        "Total_mod": R,
        "Total_phase": theta_R,
        "solution": sol,
    }

    return ret_dict


def get_max(res):
    modA = np.array([x["module_A"] for x in res])
    phiA = np.array([x["phase_A"] for x in res])

    modB = np.array([x["module_B"] for x in res])
    phiB = np.array([x["phase_B"] for x in res])

    XA = modA * np.cos(phiA)
    YA = modA * np.sin(phiA)

    XB = modB * np.cos(phiB)
    YB = modB * np.sin(phiB)

    Amax = np.argmax(modA)
    Bmax = np.argmax(modB)

    XAmax = XA[Amax]
    YAmax = YA[Amax]
    XBmax = XB[Bmax]
    YBmax = YB[Bmax]

    print(f"XAmax = {XAmax}, YAmax = {YAmax}")
    print(f"XBmax = {XBmax}, YBmax = {YBmax}")

    return 1, 1
