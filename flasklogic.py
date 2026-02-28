import numpy as np

# Constants
c = 299792.458  # Speed of light in km/s
T_CMB0 = 2.72548  # CMB temperature today in Kelvin
sec_per_Gyr = 3.15576e16  # seconds in a Gyr
Mpc_in_m = 3.085677581e22  # meters in a megaparsec

# E(z)
def E(z, Omega_m, Omega_lambda):
    Omega_r = 8.4e-5
    Omega_k = 1.0 - Omega_m - Omega_lambda - Omega_r
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_k * (1 + z)**2 + Omega_lambda)

def H(z, H0, Omega_m, Omega_lambda):
    return H0 * E(z, Omega_m, Omega_lambda)  # in km/s/Mpc

# Simpson's 3/8 rule
def simpsons_3_8(f, a, b, n, *args):
    if n % 3 != 0:
        raise ValueError("n must be a multiple of 3")
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    fx = f(x, *args)
    result = fx[0] + fx[-1]
    result += 3 * np.sum(fx[1:-1][(np.arange(1, n) % 3 != 0)])
    result += 2 * np.sum(fx[3:-3:3])
    return 3 * h * result / 8

# Integrands
def integrand_lookback(z, H0, Omega_m, Omega_lambda):
    return 1.0 / ((1 + z) * H(z, H0, Omega_m, Omega_lambda))

def integrand_comoving(z, H0, Omega_m, Omega_lambda):
    return c / H(z, H0, Omega_m, Omega_lambda)

# Main calculator
def calculate_cosmology(H0, Omega_m, Omega_lambda, z_input, n_points=300000):
    if n_points % 3 != 0:
        n_points += 3 - (n_points % 3)

    # Convert H0 to 1/s
    H0_SI = H0 * 1000 / Mpc_in_m  # km/s/Mpc → 1/s

    # Present age of the universe
    z_max = 15000
    integral_age_today = simpsons_3_8(integrand_lookback, 0, z_max, n_points, H0_SI, Omega_m, Omega_lambda)
    age_today_sec = integral_age_today 
    age_today_Gyr = age_today_sec / sec_per_Gyr

    # Lookback time
    integral_lookback = simpsons_3_8(integrand_lookback, 0, z_input, n_points, H0_SI, Omega_m, Omega_lambda)
    lookback_sec = integral_lookback
    lookback_Gyr = lookback_sec / sec_per_Gyr

    # Age at redshift z
    age_at_z_Gyr = age_today_Gyr - lookback_Gyr

    # Comoving radial distance (in Mpc)
    D_C = simpsons_3_8(integrand_comoving, 0, z_input, n_points, H0, Omega_m, Omega_lambda)

    # Angular diameter distance
    D_A = D_C / (1 + z_input)

    # Luminosity distance
    D_L = D_C * (1 + z_input)

    # Angular scale in kpc/arcsec
    angular_scale_kpc_per_arcsec = D_A * 1e3 / 206265

    # CMB temperature
    T_CMB = T_CMB0 * (1 + z_input)

    return {
        "1. Present age of the universe (in Gyr)": age_today_Gyr,
        "2. Age of the universe at given z (in Gyr)": age_at_z_Gyr,
        "3. Lookback time (in Gyr)": lookback_Gyr,
        "4. Comoving Distance (in Mpc)": D_C,
        "5. Angular Diameter Distance (in Mpc)": D_A,
        "6. Luminosity Distance (in Mpc)": D_L,
        "7. Angular Scale (in Kpc/arcsec)": angular_scale_kpc_per_arcsec,
        "8. CMB Temperature at given z (in Kelvins)": T_CMB,
    }

# Example test
if __name__ == "__main__":
    z = 2.0
    H0 = 69.6
    Omega_m = 0.286
    Omega_lambda = 0.714
    result = calculate_cosmology(H0, Omega_m, Omega_lambda, z)
    for k, v in result.items():
        print(f"{k}: {v:.4f}")
