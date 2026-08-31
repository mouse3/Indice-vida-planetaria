from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np


# -----------------------------------------------------------------------------
# 1. CONSTANTES FÍSICAS UNIVERSALES (SI)
# -----------------------------------------------------------------------------
H: float = 6.62607015e-34       # Constante de Planck [J*s]
C: float = 2.99792458e8         # Velocidad de la luz [m/s]
K_B: float = 1.380649e-23       # Constante de Boltzmann [J/K]
EV_TO_J: float = 1.602176634e-19 # Conversión de eV a Joules
R_GAS: float = 8.314462618      # Constante universal de los gases [J/(mol*K)]
SIGMA: float = 5.670374419e-8   # Constante de Stefan-Boltzmann [W/(m^2 * K^4)]


# -----------------------------------------------------------------------------
# 2. ESTROCTURAS DE DATOS Y ENLACES QUÍMICOS
# -----------------------------------------------------------------------------
@dataclass
class ChemicalBond:
    name: str
    context: str
    energy_ev: float
    energy_kj_mol: float
    lambda_nm: float


# Tabla de enlaces CHONPS y Silicio
BONDS_CHONPS: Dict[str, ChemicalBond] = {
    "OH_alcohol": ChemicalBond("O-H", "Alcoholes/Ácidos", 4.79, 463.0, 258.9),
    "OH_water":   ChemicalBond("O-H", "Agua simple",      4.81, 464.0, 257.8),
    "CC_sugar":   ChemicalBond("C-C", "Glúcidos",         3.61, 348.0, 343.5),
    "CC_single":  ChemicalBond("C-C", "Simple",           3.61, 348.0, 343.5),
    "CN_amine":   ChemicalBond("C-N", "Aminas",           3.16, 305.0, 392.4),
    "CN_single":  ChemicalBond("C-N", "Simple",           3.16, 306.0, 392.4),
    "PO_nucleo":  ChemicalBond("P-O", "Nucleótidos",      3.64, 351.0, 340.7),
    "PO_single":  ChemicalBond("P-O", "Simple",           3.64, 351.0, 340.7),
}

BONDS_SILICON: Dict[str, ChemicalBond] = {
    "SiSi_single": ChemicalBond("Si-Si", "Esqueleto silicio", 2.34, 226.0, 529.9),
    "SiO_single":  ChemicalBond("Si-O",  "Siloxanos",          4.83, 466.0, 256.7),
    "SiC_single":  ChemicalBond("Si-C",  "Organosilicones",    3.18, 307.0, 389.9),
}


# -----------------------------------------------------------------------------
# 3. RADIACIÓN Y ESPECTRO DE CUERPO NEGRO
# -----------------------------------------------------------------------------
def planck_spectral_radiance(lam_meters: np.ndarray, T_star: float) -> np.ndarray:
    """Calcula la irradiancia espectral B_lambda(lambda, T) en W/(m^2 * sr * m)."""
    exponent = (H * C) / (lam_meters * K_B * T_star)
    return (2.0 * H * C**2 / (lam_meters**5)) / (np.exp(exponent) - 1.0)


def photon_energy_ev(lam_meters: np.ndarray) -> np.ndarray:
    """Convierte longitud de onda en metros a energía de fotón en eV."""
    return (H * C / lam_meters) / EV_TO_J


def calculate_surface_flux(T_star: float, R_star: float, distance_au: float, 
                           lam_range: Tuple[float, float]) -> float:
    """Calcula el flujo de radiación espectral integrado [W/m^2] a una distancia planetaria d."""
    AU_in_meters = 1.496e11
    d_meters = distance_au * AU_in_meters
    
    lam = np.linspace(lam_range[0], lam_range[1], 1000)
    b_lam = planck_spectral_radiance(lam, T_star)
    
    # Integración trapezoidal del espectro estelar
    integrated_radiance = np.trapezoid(b_lam, lam) # [W / (m^2 sr)]
    flux_at_star = np.pi * integrated_radiance # [W / m^2]
    
    # Ley de la inversa del cuadrado
    flux_at_planet = flux_at_star * (R_star / d_meters)**2
    return flux_at_planet


# -----------------------------------------------------------------------------
# 4. VENTANA FOTOQUÍMICA Y PROBABILIDAD DE ACTIVACIÓN
# -----------------------------------------------------------------------------
def single_bond_probability(E_ev: np.ndarray, 
                            E_b_ev: float, 
                            E_exc_ev: float = 1.8, 
                            beta1: float = 10.0, 
                            beta2: float = 10.0) -> np.ndarray:
    """Calcula W(E, E_b) para un enlace específico."""
    p_activation = 1.0 / (1.0 + np.exp(-beta1 * (E_ev - E_exc_ev)))
    p_survival = 1.0 / (1.0 + np.exp(beta2 * (E_ev - E_b_ev)))
    return p_activation * p_survival


def total_photochemical_window(E_ev: np.ndarray, 
                               bonds: Dict[str, ChemicalBond], 
                               E_exc_ev: float = 1.8, 
                               beta1: float = 10.0, 
                               beta2: float = 10.0) -> np.ndarray:
    """
    Calcula la ventana fotoquímica total W_T(E) mediante la unión probabilística:
    W_T(E) = 1 - PROD_i (1 - W_i(E, E_b_i)) en el intervalo [0, 1].
    """
    prod_term = np.ones_like(E_ev, dtype=float)
    for bond in bonds.values():
        w_i = single_bond_probability(E_ev, bond.energy_ev, E_exc_ev, beta1, beta2)
        prod_term *= (1.0 - w_i)
    
    return 1.0 - prod_term


def photoelectric_kinetic_energy(E_ev: np.ndarray, work_function_ev: float) -> np.ndarray:
    """Calcula E_c_max = h*nu - W_0 (en eV), limitando valores negativos a 0."""
    return np.maximum(0.0, E_ev - work_function_ev)


# -----------------------------------------------------------------------------
# 5. TERMODINÁMICA Y ENTROPÍA DEL FLUJO DE FOTONES
# -----------------------------------------------------------------------------
def arrhenius_rate_constant(A: float, m: float, E_a_J: float, T_planet: float) -> float:
    """Calcula la constante cinética de reacción: k = A * T^m * exp(-E_a / (R * T))."""
    return A * (T_planet**m) * np.exp(-E_a_J / (R_GAS * T_planet))


def negentropy_production_rate(E_absorbed_joules: float, 
                               T_planet: float, 
                               T_star: float) -> float:
    """
    Calcula la tasa de producción de negentropía disponible:
    Delta S = (4/3) * E * (1 / T_planet - 1 / T_star)  [J / (K * m^2 * s)]
    """
    return (4.0 / 3.0) * E_absorbed_joules * ((1.0 / T_planet) - (1.0 / T_star))


# -----------------------------------------------------------------------------
# 6. SIMULACIÓN INTEGRAL DEL MODELO
# -----------------------------------------------------------------------------
class AstrobiologicalHabitabilityModel:
    def __init__(self, T_star: float, R_star: float, distance_au: float, T_planet: float):
        self.T_star = T_star
        self.R_star = R_star
        self.distance_au = distance_au
        self.T_planet = T_planet

    def evaluate(self, bonds: Dict[str, ChemicalBond], E_exc_ev: float = 1.8) -> dict:
        # Rango espectral (100 nm a 2000 nm)
        lam_grid = np.linspace(100e-9, 2000e-9, 1000)
        E_grid_ev = photon_energy_ev(lam_grid)
        
        # Radiación de cuerpo negro e irradiancia espectral
        flux_total = calculate_surface_flux(self.T_star, self.R_star, self.distance_au, (100e-9, 2000e-9))
        
        # Ventana fotoquímica integrada
        w_t = total_photochemical_window(E_grid_ev, bonds, E_exc_ev=E_exc_ev)
        photochemical_habitability_index = float(np.trapezoid(w_t, lam_grid) / (lam_grid[-1] - lam_grid[0]))
        
        # Balance Negentrópico
        delta_S = negentropy_production_rate(flux_total, self.T_planet, self.T_star)

        return {
            "flux_total_W_m2": flux_total,
            "photochemical_index": photochemical_habitability_index,
            "negentropy_production_J_K_m2_s": delta_S
        }


# Ejemplo de uso
if __name__ == "__main__":
    # Caso 1: Tierra orbitando al Sol
    sol_tierra = AstrobiologicalHabitabilityModel(
        T_star=3547.0, 
        R_star=6.9634e8*0.41, 
        distance_au=0.159, 
        T_planet=265.0
    )
    
    res_chonps = sol_tierra.evaluate(BONDS_CHONPS, E_exc_ev=1.8)
    res_silicon = sol_tierra.evaluate(BONDS_SILICON, E_exc_ev=1.5)

    print("--- RESULTADOS SIMULACIÓN CHONPS ---")
    print(f"Flujo Total Absorbido: {res_chonps['flux_total_W_m2']:.2f} W/m^2")
    print(f"Índice de Ventana Fotoquímica W_T: {res_chonps['photochemical_index']:.4f}")
    print(f"Producción de Negentropía (Delta S): {res_chonps['negentropy_production_J_K_m2_s']:.4f} J/(K*m^2*s)")

    print("\n--- RESULTADOS SIMULACIÓN SILICIO ---")
    print(f"Índice de Ventana Fotoquímica W_T (Silicio): {res_silicon['photochemical_index']:.4f}")