import pathlib
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict
from scipy.special import expit

PLOT_FOLDER = pathlib.Path(__file__).resolve().parent.parent / "plots"


def plot_piston():
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    t = np.linspace(0, 1, 100) * 2 * np.pi
    y = expit((t - np.pi) * 20)

    ax.plot(t, y)

    ax.set_xlabel("time [s]")
    ax.set_ylabel("Step Response")
    ax.grid(True)

    fig.savefig(PLOT_FOLDER / "step_response.png")


def plot_results(res: List[Dict[str, np.float32]], total_time):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    total_mod = np.array([x["Total_mod"] for x in res])
    total_phase = np.array([x["Total_phase"] for x in res])*180/np.pi

    ax.plot(total_time, total_mod, label="Resultante Total [kN]")
    ax.plot(total_time, total_phase, label="Fase Total [Deg]")

    ax.set_xlabel("Giro de la manivela [rad]")
    # ax.set_ylabel("Modulo [kN] / Phase [Deg]")
    ax.grid(True)
    ax.legend()

    fig.savefig(PLOT_FOLDER / "total_mod.png")

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    XAb = np.array([x["solution"][0] for x in res])
    YAb = np.array([x["solution"][1] for x in res])

    ax.scatter(XAb, YAb)

    ax.set_xlabel("XAb [kN]")
    ax.set_ylabel("YAb [kN]")
    ax.grid(True)

    fig.savefig(PLOT_FOLDER / "reactionsA.png")

    fig, ax = plt.subplots(2, 1, figsize=(12, 8))

    modA = np.array([x["module_A"] for x in res])
    pA = np.array([x["phase_A"] for x in res])*180/np.pi
    modB = np.array([x["module_B"] for x in res])
    pB = np.array([x["phase_B"] for x in res])*180/np.pi

    ax[0].plot(total_time, modA, label="Modulo [kN]")
    ax[1].plot(total_time, modB, label="Modulo [kN]")
    ax[0].plot(total_time, pA, label="Fase A [Deg]")
    ax[1].plot(total_time, pB, label="Fase B [Deg]")

    ax[0].grid(True)
    ax[0].legend()
    ax[1].grid(True)
    ax[1].legend()

    ax[0].set_xlabel("Giro de la manivela [rad]")
    ax[1].set_xlabel("Giro de la manivela [rad]")
    ax[0].set_ylabel("Modulo Fuerza Punto A [kN]")
    ax[1].set_ylabel("Modulo Fuerza Punto B [kN]")

    fig.savefig(PLOT_FOLDER / "modules.png")