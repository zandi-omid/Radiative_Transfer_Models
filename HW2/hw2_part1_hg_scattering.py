"""HW2 baseline: HG scattering without absorption."""

from models.monte_carlo import MonteCarloTransport


def main():
    n_photons = 10_000
    tau = 4.0
    mu0 = -0.7
    g = 0.75
    omega0 = 1.0

    model = MonteCarloTransport(
        tau=tau,
        mu0=mu0,
        g=g,
        omega0=omega0,
        workers=24,
    )
    result = model.run(n_photons, seed=42)

    print("HW2 - Baseline HG Case")
    print("----------------------")
    print(f"Photons             : {n_photons:,}")
    print(f"Optical depth       : {tau}")
    print(f"Initial mu          : {mu0}")
    print(f"Asymmetry factor g  : {g}")
    print(f"Single scat. albedo : {omega0}")
    print()
    print(f"Reflectivity        : {result.reflectivity:.6f}")
    print(f"Transmissivity      : {result.transmissivity:.6f}")
    print(f"Absorptivity        : {result.absorptivity:.6f}")
    print(
        "R + T + A           : "
        f"{result.reflectivity + result.transmissivity + result.absorptivity:.6f}"
    )
    print(f"Mean scatterings    : {result.mean_scatterings:.3f}")
    print(f"Elapsed time        : {result.elapsed_seconds:.3f} s")


if __name__ == "__main__":
    main()
