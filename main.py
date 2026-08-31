import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple

# -----------------------------------------------------------------------------
# 1. CONSTANTES FÍSICAS UNIVERSALES (SI)
# -----------------------------------------------------------------------------
H: float = 6.62607015e-34       # Constante de Planck [J*s]
C: float = 2.99792458e8         # Velocidad de la luz [m/s]
K_B: float = 1.380649e-23       # Constante de Boltzmann [J/K]
EV_TO_J: float = 1.602176634e-19 # Conversión de eV a Joules


# -----------------------------------------------------------------------------
# 2. ESTRUCTURA DE DATOS Y TABLAS DE ENLACES QUÍMICOS
# -----------------------------------------------------------------------------
@dataclass
class ChemicalBond:
    name: str
    context: str
    energy_ev: float
    energy_kj_mol: float
    lambda_nm: float


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
# 3. FUNCIONES AUXILIARES DE RADIACIÓN Y PROBABILIDAD FOTOQUÍMICA
# -----------------------------------------------------------------------------
def planck_spectral_radiance(lam_meters: np.ndarray, T_star: float) -> np.ndarray:
    """Calcula la irradiancia espectral B_lambda(lambda, T) en W/(m^2 * sr * m)."""
    exponent = (H * C) / (lam_meters * K_B * T_star)
    return (2.0 * H * C**2 / (lam_meters**5)) / (np.exp(exponent) - 1.0)


def photon_energy_ev(lam_meters: np.ndarray) -> np.ndarray:
    """Convierte longitud de onda en metros a energía de fotón en eV."""
    return (H * C / lam_meters) / EV_TO_J


def single_bond_probability(E_ev: np.ndarray, E_b_ev: float, E_exc_ev: float = 1.8, 
                            beta1: float = 10.0, beta2: float = 10.0) -> np.ndarray:
    """Calcula la probabilidad W_i(E, E_b) para un enlace individual."""
    p_activation = 1.0 / (1.0 + np.exp(-beta1 * (E_ev - E_exc_ev)))
    p_survival = 1.0 / (1.0 + np.exp(beta2 * (E_ev - E_b_ev)))
    return p_activation * p_survival


def total_photochemical_window(E_ev: np.ndarray, bonds: Dict[str, ChemicalBond], 
                               E_exc_ev: float = 1.8, beta1: float = 10.0, beta2: float = 10.0) -> np.ndarray:
    """Unión probabilística de ventanas individuales: W_T(E) = 1 - PROD(1 - W_i)."""
    prod_term = np.ones_like(E_ev, dtype=float)
    for bond in bonds.values():
        w_i = single_bond_probability(E_ev, bond.energy_ev, E_exc_ev, beta1, beta2)
        prod_term *= (1.0 - w_i)
    return 1.0 - prod_term


# -----------------------------------------------------------------------------
# 4. CLASE PRINCIPAL DEL MODELO DE HABITABILIDAD ASTROBIOLÓGICA
# -----------------------------------------------------------------------------
class AstrobiologicalHabitabilityModel:
    def __init__(self, T_star: float, R_star: float, distance_au: float, 
                 R_planet: float, T_planet: float, albedo: float = 0.3):
        """
        :param T_star: Temperatura efectiva de la estrella [K]
        :param R_star: Radio de la estrella [m]
        :param distance_au: Distancia orbital [UA]
        :param R_planet: Radio del planeta [m]
        :param T_planet: Temperatura media superficial del planeta [K]
        :param albedo: Albedo de Bond planetario [0 a 1]
        """
        self.T_star = T_star
        self.R_star = R_star
        self.distance_m = distance_au * 1.496e11
        self.R_planet = R_planet
        self.T_planet = T_planet
        self.albedo = albedo

        # Geometría del planeta
        self.A_cross_section = np.pi * (R_planet ** 2)  # Disco captador de radiación [m^2]
        self.A_surface = 4.0 * np.pi * (R_planet ** 2)   # Esfera planetaria total [m^2]

    def evaluate(self, bonds: Dict[str, ChemicalBond], E_exc_ev: float = 1.8, 
                 lam_range: Tuple[float, float] = (100e-9, 2000e-9), num_points: int = 1000) -> dict:
        
        lam_grid = np.linspace(lam_range[0], lam_range[1], num_points)
        E_grid_ev = photon_energy_ev(lam_grid)

        # 1. Radiación espectral de cuerpo negro
        b_lam = planck_spectral_radiance(lam_grid, self.T_star)
        
        # 2. Irradiancia en la superficie estelar y en la órbita (E/A incidente)
        flux_star_surface = np.pi * np.trapz(b_lam, lam_grid)
        flux_incidente = flux_star_surface * (self.R_star / self.distance_m)**2

        # 3. Potencia Interceptada y Absorbida considerando R_planeta y Albedo
        power_intercepted = flux_incidente * self.A_cross_section
        power_absorbed = power_intercepted * (1.0 - self.albedo)
        flux_surface_avg = power_absorbed / self.A_surface

        # 4. Ventana fotoquímica ponderada espectralmente por la emisión de la estrella
        w_t_lam = total_photochemical_window(E_grid_ev, bonds, E_exc_ev=E_exc_ev)
        
        total_spectral_energy = np.trapz(b_lam, lam_grid)
        useful_spectral_energy = np.trapz(w_t_lam * b_lam, lam_grid)
        
        # Índice de compatibilidad fotoquímica normalizado [0, 1]
        photochemical_index = float(useful_spectral_energy / total_spectral_energy)

        # 5. Tasa global de generación de negentropía disipativa
        negentropy_rate_global = (4.0 / 3.0) * power_absorbed * ((1.0 / self.T_planet) - (1.0 / self.T_star))

        return {
            "flux_incidente_W_m2": flux_incidente,
            "power_intercepted_GW": power_intercepted / 1e9,
            "power_absorbed_GW": power_absorbed / 1e9,
            "flux_surface_avg_W_m2": flux_surface_avg,
            "photochemical_index": photochemical_index,
            "negentropy_rate_global_MW_K": negentropy_rate_global / 1e6
        }


# -----------------------------------------------------------------------------
# 5. EJECUCIÓN Y COMPARACIÓN DE ESCENARIOS (TIERRA VS K2-18b)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    R_SUN = 6.9634e8
    R_EARTH = 6.371e6

    # Configuración Sistema Solar - Tierra
    tierra_model = AstrobiologicalHabitabilityModel(
        T_star=5778.0,
        R_star=R_SUN,
        distance_au=1.0,
        R_planet=R_EARTH,
        T_planet=288.0,
        albedo=0.3
    )

    # Configuración K2-18 - K2-18b
    k218b_model = AstrobiologicalHabitabilityModel(
        T_star=3547.0,
        R_star=R_SUN * 0.41,
        distance_au=0.159,
        R_planet=R_EARTH * 2.61,
        T_planet=265.0,
        albedo=0.3
    )

    # Evaluación CHONPS vs Silicio
    tierra_chonps = tierra_model.evaluate(BONDS_CHONPS, E_exc_ev=1.8)
    tierra_silicon = tierra_model.evaluate(BONDS_SILICON, E_exc_ev=1.5)

    k218b_chonps = k218b_model.evaluate(BONDS_CHONPS, E_exc_ev=1.8)
    k218b_silicon = k218b_model.evaluate(BONDS_SILICON, E_exc_ev=1.5)

    print("==================================================")
    print("               RESULTADOS: TIERRA                 ")
    print("==================================================")
    print(f"Flujo Incidente (E/A):               {tierra_chonps['flux_incidente_W_m2']:.2f} W/m^2")
    print(f"Flujo Promedio Superficie (Abs):      {tierra_chonps['flux_surface_avg_W_m2']:.2f} W/m^2")
    print(f"Índice W_T Ponderado (CHONPS):        {tierra_chonps['photochemical_index']:.4f}")
    print(f"Índice W_T Ponderado (Silicio):       {tierra_silicon['photochemical_index']:.4f}")
    print(f"Generación Negentropía Global:       {tierra_chonps['negentropy_rate_global_MW_K']:,.2f} MW/K\n")

    print("==================================================")
    print("              RESULTADOS: K2-18b                  ")
    print("==================================================")
    print(f"Flujo Incidente (E/A):               {k218b_chonps['flux_incidente_W_m2']:.2f} W/m^2")
    print(f"Flujo Promedio Superficie (Abs):      {k218b_chonps['flux_surface_avg_W_m2']:.2f} W/m^2")
    print(f"Índice W_T Ponderado (CHONPS):        {k218b_chonps['photochemical_index']:.4f}")
    print(f"Índice W_T Ponderado (Silicio):       {k218b_silicon['photochemical_index']:.4f}")
    print(f"Generación Negentropía Global:       {k218b_chonps['negentropy_rate_global_MW_K']:,.2f} MW/K\n")