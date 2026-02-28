import streamlit as st
import numpy as np

# Constants
c = 299792.458  # Speed of light in km/s
T_CMB0 = 2.72548  # CMB temperature today in Kelvin
sec_per_Gyr = 3.15576e16  # seconds in a Gyr
Mpc_in_m = 3.085677581e22  # meters in a megaparsec

st.set_page_config(page_title="Sourya's Python Cosmology Calculator", layout="wide")
st.title("🌌 Cosmology Calculator")

st.markdown("""
This calculator is inspired by [Ned Wright's Cosmology Calculator](http://www.astro.ucla.edu/~wright/CosmoCalc.html).
Use it to compute cosmological distances, time scales, and CMB temperature at a given redshift.

Modify the default cosmological parameters if needed and input the desired redshift.
""")

# Sidebar for input
with st.sidebar:
    st.header("🔧 Cosmological Parameters")
    H0 = st.number_input("Hubble Constant H₀ (km/s/Mpc)", value=69.6, min_value=0.0, step=0.1)
    Omega_m = st.number_input("Matter Density Ωₘ₀", value=0.286, min_value=0.0, max_value=2.0, step=0.001, format="%.3f")
    Omega_lambda = st.number_input("Dark Energy Density Ω_Λ₀", value=0.714, min_value=0.0, max_value=2.0, step=0.001, format="%.3f")
    Omega_k = 1.0 - Omega_m - Omega_lambda
    st.write(f"Curvature Density Ωₖ = {Omega_k:.5f} (auto-computed)")

    # Universe geometry interpretation
    if round(Omega_k, 2) == 0.00:
        st.success("Universe is **Flat** (Ωₖ ≈ 0)")
    elif round(Omega_k, 2) < 0.00:
        st.warning("Universe is **Closed** (Ωₖ < 0)")
    else:
        st.info("Universe is **Open** (Ωₖ > 0)")

    z = st.number_input("Redshift z", min_value=0.0, value=1.0, step=0.01)

# Define functions
def E(z, Omega_m, Omega_lambda):
    Omega_r = 8.4e-5
    Omega_k = 1.0 - Omega_m - Omega_lambda
    return np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_k * (1 + z)**2 + Omega_lambda)

def H(z, H0, Omega_m, Omega_lambda):
    return H0 * E(z, Omega_m, Omega_lambda)

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

def integrand_lookback(z, H0_SI, Omega_m, Omega_lambda):
    return 1.0 / ((1 + z) * H(z, H0_SI, Omega_m, Omega_lambda))

def integrand_comoving(z, H0, Omega_m, Omega_lambda):
    return c / H(z, H0, Omega_m, Omega_lambda)

def calculate_cosmology(H0, Omega_m, Omega_lambda, z_input, n_points=300000):
    if n_points % 3 != 0:
        n_points += 3 - (n_points % 3)

    H0_SI = H0 * 1000 / Mpc_in_m  # Convert to 1/s

    z_max = 15000
    age_today_sec = simpsons_3_8(integrand_lookback, 0, z_max, n_points, H0_SI, Omega_m, Omega_lambda)
    age_today_Gyr = age_today_sec / sec_per_Gyr

    lookback_sec = simpsons_3_8(integrand_lookback, 0, z_input, n_points, H0_SI, Omega_m, Omega_lambda)
    lookback_Gyr = lookback_sec / sec_per_Gyr
    age_at_z_Gyr = age_today_Gyr - lookback_Gyr

    D_C = simpsons_3_8(integrand_comoving, 0, z_input, n_points, H0, Omega_m, Omega_lambda)
    D_A = D_C / (1 + z_input)
    D_L = D_C * (1 + z_input)
    angular_scale_kpc_per_arcsec = D_A * 1e3 / 206265

    if z_input <= 1249:
        T_CMB = T_CMB0 * (1 + z_input)
    else:
        T_CMB = None  # No CMB before recombination

    return {
        "age_today": age_today_Gyr,
        "age_at_z": age_at_z_Gyr,
        "lookback_time": lookback_Gyr,
        "comoving_distance_Mpc": D_C,
        "angular_diameter_distance_Mpc": D_A,
        "luminosity_distance_Mpc": D_L,
        "angular_scale_kpc_per_arcsec": angular_scale_kpc_per_arcsec,
        "CMB_temperature_K": T_CMB,
    }

# Compute and display
if st.button("📡 Compute Cosmology"):
    results = calculate_cosmology(H0, Omega_m, Omega_lambda, z)

    st.subheader(f"Results at z = {z}")
    col1, col2 = st.columns(2)

    # Assign time unit and factor individually
    age_today_unit = "Gyr" if results["age_today"] >= 1 else "Myr"
    age_today_factor = 1.0 if results["age_today"] >= 1 else 1000.0

    age_at_z_unit = "Gyr" if results["age_at_z"] >= 1 else "Myr"
    age_at_z_factor = 1.0 if results["age_at_z"] >= 1 else 1000.0

    lookback_unit = "Gyr" if results["lookback_time"] >= 1 else "Myr"
    lookback_factor = 1.0 if results["lookback_time"] >= 1 else 1000.0

    with col1:
        st.metric(f"🕰 Present Age of Universe", f"{results['age_today'] * age_today_factor:.5f} {age_today_unit}")
        st.metric(f"🕰 Age of Universe at z", f"{results['age_at_z'] * age_at_z_factor:.5f} {age_at_z_unit}")
        st.metric(f"⏳ Lookback Time", f"{results['lookback_time'] * lookback_factor:.5f} {lookback_unit}")
        st.metric("📏 Comoving Radial Distance", f"{results['comoving_distance_Mpc']:.5f} Mpc")

    with col2:
        st.metric("📐 Angular Diameter Distance", f"{results['angular_diameter_distance_Mpc']:.5f} Mpc")
        st.metric("🧭 Angular Scale", f"{results['angular_scale_kpc_per_arcsec']:.5f} kpc/arcsec")
        st.metric("💡 Luminosity Distance", f"{results['luminosity_distance_Mpc']:.5f} Mpc")

        if results["CMB_temperature_K"] is not None:
            st.metric("🌡 CMB Temperature at z", f"{results['CMB_temperature_K']:.5f} K")
        else:
            st.metric("🌡 CMB Temperature at z", "N/A (before recombination)")

    st.markdown("---")
    st.caption("""
    📎 **Footnotes**:
    - Speed of light used: c = 299,792.458 km/s 
    - 1 Gyr = 1,000,000,000 years 
    - 1 Myr = 1,000,000 years 
    - 1 Mpc = 1,000,000 parsecs 
    - 1 parsec = 3.08568 × 10¹⁶ m 
    - Default values: Ωₘ = 0.286, Ω_Λ = 0.714, Ωᵣ = 8.4 × 10⁻⁵, Ωₖ = 0
    """)

else:
    st.info("Enter a redshift and press **Compute Cosmology**.")
