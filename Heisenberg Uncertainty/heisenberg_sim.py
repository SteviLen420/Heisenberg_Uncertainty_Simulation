# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Stefan Len
#
# ===================================================================================
# Heisenberg_Uncertainty_Minimal_Simulation.py
# ===================================================================================
# Author: Stefan Len
# Overview:
#   Minimal numerical demonstration of the Heisenberg Uncertainty Principle (Colab-ready):
#   - Auto-mounts Google Drive (if in Colab)
#   - Installs missing Python packages quietly
#   - Generates Gaussian wave packets in position space
#   - Computes uncertainties Δx and Δp via FFT
#   - Verifies the uncertainty relation (Δx·Δp ≥ ħ/2)
#   - Optional free-particle time evolution (wave packet spreading)
#   - Saves ALL outputs to Drive under MyDrive/Heisenberg_Sims/run_YYYYMMDD_HHMMSS/
# ===================================================================================

import sys, subprocess, importlib, os, datetime

# ---------- Helper: silent pip install if a module is missing ----------
def ensure_package(mod_name, pip_name=None):
    try:
        importlib.import_module(mod_name)
    except ModuleNotFoundError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", pip_name or mod_name])

# Ensure core deps (usually present in Colab, but make it robust elsewhere)
ensure_package("numpy")
ensure_package("matplotlib")

# Now safe to import
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from numpy.fft import fft, fftshift, ifft
import csv

# ---------- Detect Colab & mount Drive ----------
def in_colab():
    try:
        import google.colab  # noqa: F401
        return True
    except Exception:
        return False

DRIVE_MOUNTED = False
BASE_DIR = os.path.abspath("./heisenberg_outputs_local")

if in_colab():
    try:
        from google.colab import drive
        drive.mount("/content/drive", force_remount=True)
        DRIVE_MOUNTED = True
        BASE_DIR = "/content/drive/MyDrive/Heisenberg_Sims"
    except Exception as e:
        print("[WARN] Could not mount Google Drive. Falling back to local folder.", e)
        DRIVE_MOUNTED = False
        BASE_DIR = os.path.abspath("./heisenberg_outputs_local")
else:
    # Not Colab (e.g. local Jupyter/Python) → use local folder
    BASE_DIR = os.path.abspath("./heisenberg_outputs_local")

# ---------- Create run directory with timestamp ----------
RUN_ID = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
outdir = os.path.join(BASE_DIR, f"run_{RUN_ID}")
os.makedirs(outdir, exist_ok=True)

# ---------- Save a small run metadata file ----------
with open(os.path.join(outdir, "run_info.txt"), "w", encoding="utf-8") as fmeta:
    fmeta.write(
        "Heisenberg Uncertainty Minimal Simulation (Colab-ready)\n"
        f"RUN_ID: {RUN_ID}\n"
        f"BASE_DIR: {BASE_DIR}\n"
        f"OUTDIR: {outdir}\n"
        f"Python: {sys.version.split()[0]}\n"
    )

# ===== Fundamental constants & grids =====
hbar = 1.0        # unit system: ℏ = 1
m    = 1.0        # particle mass
N    = 2**14      # number of grid points (large enough but still efficient)
L    = 200.0      # spatial domain: x ∈ [-L/2, L/2)
dx   = L / N
x    = (np.arange(N) - N//2) * dx

# Fourier grid (momentum space)
dk = 2*np.pi / L
k  = (np.arange(N) - N//2) * dk
p  = hbar * k

def normalize(psi, dx):
    """Normalize a wavefunction in x-space."""
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
    return psi / norm

def gaussian_packet(x, x0=0.0, p0=0.0, sigma_x=2.0):
    """Gaussian wave packet centered at x0 with mean momentum p0."""
    A = (1.0/(2*np.pi*sigma_x**2))**0.25
    return A * np.exp(- (x-x0)**2 / (4*sigma_x**2) + 1j * p0 * (x-x0)/hbar)

def expectation_and_sigma(x, density, dx):
    """Compute expectation value and standard deviation for a 1D density over x."""
    Ex  = np.sum(x * density) * dx
    Ex2 = np.sum(x**2 * density) * dx
    var = max(Ex2 - Ex**2, 0.0)
    return Ex, np.sqrt(var)

def sigma_p_from_fft(psi, dx):
    """Compute momentum expectation and spread using Fourier transform."""
    Psi_k  = fftshift(fft(np.fft.ifftshift(psi))) * dx/np.sqrt(2*np.pi)  # unitary conv.
    dens_p = np.abs(Psi_k)**2
    dp     = p[1] - p[0]
    dens_p /= np.sum(dens_p) * dp  # normalize in p-space
    Ep, sigp = expectation_and_sigma(p, dens_p, dp)
    return Ep, sigp, dens_p

def main():
    # ===== 1) Single Gaussian example (plots) =====
    sigma_x0 = 3.0
    p0       = 0.0
    psi      = gaussian_packet(x, x0=0.0, p0=p0, sigma_x=sigma_x0)
    psi      = normalize(psi, dx)

    dens_x = np.abs(psi)**2
    Ex, sigx = expectation_and_sigma(x, dens_x, dx)
    Ep, sigp, dens_p = sigma_p_from_fft(psi, dx)

    # Position-space density
    plt.figure(figsize=(7,4.5))
    plt.plot(x, dens_x)
    plt.xlabel("x")
    plt.ylabel(r"$|\psi(x)|^2$")
    plt.title(fr"Position density, $\sigma_x \approx {sigx:.3f}$")
    plt.xlim(-30, 30)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"position_density_sigma{sigma_x0:.1f}.png"), dpi=300)

    # Momentum-space density
    plt.figure(figsize=(7,4.5))
    plt.plot(p, dens_p)
    plt.xlabel("p")
    plt.ylabel(r"$|\Psi(p)|^2$")
    plt.title(fr"Momentum density, $\sigma_p \approx {sigp:.3f}$")
    plt.xlim(-3, 3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"momentum_density_sigma{sigma_x0:.1f}.png"), dpi=300)

    # Console summary
    print(f"[Single Gaussian] sigma_x ≈ {sigx:.4f}, sigma_p ≈ {sigp:.4f}, "
          f"product ≈ {sigx*sigp:.4f} (ħ/2 = {hbar/2:.4f})")

    # ===== 2) Sweep sigma_x values and measure uncertainty product =====
    sigmas = np.geomspace(0.5, 8.0, 16)
    rows = []
    for sx in sigmas:
        psi = gaussian_packet(x, sigma_x=sx)
        psi = normalize(psi, dx)
        dens_x = np.abs(psi)**2
        _, sigx_m = expectation_and_sigma(x, dens_x, dx)
        _, sigp_m, _ = sigma_p_from_fft(psi, dx)
        rows.append([sx, sigx_m, sigp_m, sigx_m*sigp_m])

    csv_path = os.path.join(outdir, "heisenberg_scan.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sigma_x_input", "sigma_x_measured", "sigma_p_measured", "product_sigma"])
        w.writerows(rows)

    # Uncertainty relation plot
    rows = np.array(rows)
    plt.figure(figsize=(7,4.5))
    plt.plot(rows[:,1], rows[:,3], "o-")
    plt.axhline(hbar/2, color="k", linestyle="--", label=r"$\hbar/2$")
    plt.xlabel(r"$\sigma_x$")
    plt.ylabel(r"$\sigma_x \sigma_p$")
    plt.title("Uncertainty product vs. position spread")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "uncertainty_product.png"), dpi=300)

    # ===== 3) (Optional) Free evolution using split-operator =====
    T  = 1.0
    dt = 0.002
    steps = int(T/dt)
    psi_t = gaussian_packet(x, sigma_x=2.0)
    psi_t = normalize(psi_t, dx)

    K = (p**2)/(2*m)
    phase_half = np.exp(-1j * K * dt/(2*hbar))

    sigx_t = []
    for _ in range(steps):
        # half step in p-space
        Psi = fftshift(fft(np.fft.ifftshift(psi_t))) * dx/np.sqrt(2*np.pi)
        Psi *= phase_half
        psi_t = fftshift(ifft(np.fft.ifftshift(Psi))) * (np.sqrt(2*np.pi)/dx)

        # (no potential in x-space)

        # second half step in p-space
        Psi = fftshift(fft(np.fft.ifftshift(psi_t))) * dx/np.sqrt(2*np.pi)
        Psi *= phase_half
        psi_t = fftshift(ifft(np.fft.ifftshift(Psi))) * (np.sqrt(2*np.pi)/dx)

        dens = np.abs(psi_t)**2
        _, sx_t = expectation_and_sigma(x, dens, dx)
        sigx_t.append(sx_t)

    # Plot free spreading of Gaussian packet
    plt.figure(figsize=(7,4.5))
    plt.plot(np.arange(steps)*dt, sigx_t)
    plt.xlabel("time")
    plt.ylabel(r"$\sigma_x(t)$")
    plt.title("Free spreading of a Gaussian wave packet")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "free_spreading_sigma_x_t.png"), dpi=300)

    print("\nOK. All figures and CSV are stored in:")
    print(outdir)
    if DRIVE_MOUNTED:
        print("Google Drive mounted: YES (MyDrive/Heisenberg_Sims)")
    else:
        print("Google Drive mounted: NO (saved locally).")

if __name__ == "__main__":
    main()
