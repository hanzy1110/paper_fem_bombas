# import fire
import numpy as np
from src.solve import solve_system, get_max
from src.plots import plot_piston, plot_results


def main():
    n_steps = 10000

    L = 330 * 1e-3  # m
    Lg = 134 * 1e-3  # m

    mb = 23514.40 * 1e-3  # kg
    mp = 23514.40 * 1e-3  # kg
    # Ib = 85570845.08/(1000**3)  # g mm**2
    Ib = 85570845.08 * 1e-9  # kgm**2
    R = 3 * 25.4 * 1e-3  # mm
    D_plunger = 114.3 * 1e-3  # mm
    p_linea = 0.0 * 1e5  # N/m**2
    # p_servicio = 443* 1e5  # N/m**2
    p_servicio = 117 * 1e5  # N/m**2

    # omega = (230 * 2 * np.pi) / 60  # rad/s
    omega = (2070 * 2 * np.pi) / 60  # rad/s
    ts = np.linspace(0, 1, n_steps) * (10 / omega)

    res = []
    for t in ts:
        alfa = (omega * t) % (2*np.pi)
        beta = np.sqrt(1 - (R / L) ** 2 * np.sin(alfa))

        r = solve_system(
            alfa, beta, omega, mb, mp, R, Ib, L, Lg, D_plunger, p_linea, p_servicio
        )
        res.append(r)

    plot_results(res, ts*omega)
    plot_piston()
    Xmax, Ymax = get_max(res)


if __name__ == "__main__":
    main()
    # fire.Fire(main)
