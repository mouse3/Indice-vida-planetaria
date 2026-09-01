import numpy as np
from scipy.optimize import fsolve
from dataclasses import dataclass
from typing import Dict, Tuple, Optional


# 1. CONSTANTES FÍSICAS UNIVERSALES (SI)
H: float = 6.62607015e-34       # Constante de Planck [J*s]
C: float = 2.99792458e8         # Velocidad de la luz [m/s]
K_B: float = 1.380649e-23       # Constante de Boltzmann [J/K]
EV_TO_J: float = 1.602176634e-19 # Conversión de eV a Joules
R_GAS: float = 8.314462618      # Constante universal de los gases [J/(mol*K)]

R_SUN = 6.9634e8                # Radio solar [m]
R_EARTH = 6.371e6               # Radio terrestre [m]
AU_M = 1.496e11                 # Unidad Astronómica [m]



# 2. ESTRUCTURA DE DATOS (ENLACES, REACCIONES Y PLANETAS)
@dataclass
class ChemicalBond:
    name: str
    context: str
    energy_ev: float
    energy_kj_mol: float
    lambda_nm: float

@dataclass
class BiokineticReaction:
    name: str
    context: str
    A_pre_exp: float
    E_a_kj_mol: float
    dH_denat_kj_mol: float
    T_m_k: float

@dataclass
class PlanetConfig:
    name: str
    T_star: float
    R_star_sun: float
    distance_au: float
    R_planet_earth: float
    T_planet: float
    albedo: float
    eccentricity: float
    rotation_period_s: float
    orbital_period_s: float
    theta_deg: float
    phi_deg: float

# Tablas de Química y Biocinética
BONDS_CHONPS: Dict[str, ChemicalBond] = {
    "OH_alcohol": ChemicalBond("O-H", "Alcoholes/Ácidos", 4.79, 463.0, 258.9),
    "OH_water":   ChemicalBond("O-H", "Agua simple",      4.81, 464.0, 257.8),
    "CC_single":  ChemicalBond("C-C", "Simple",           3.61, 348.0, 343.5),
    "CN_single":  ChemicalBond("C-N", "Simple",           3.16, 306.0, 392.4),
    "PO_single":  ChemicalBond("P-O", "Simple",           3.64, 351.0, 340.7),
}

BONDS_SILICON: Dict[str, ChemicalBond] = {
    "SiSi_single": ChemicalBond("Si-Si", "Esqueleto",    2.34, 226.0, 529.9),
    "SiO_single":  ChemicalBond("Si-O",  "Siloxanos",    4.83, 466.0, 256.7),
    "SiC_single":  ChemicalBond("Si-C",  "Organosilicio",3.18, 307.0, 389.9),
}

BIO_REACTIONS: Dict[str, BiokineticReaction] = {
    "atp_hydrolysis": BiokineticReaction("Hidrólisis ATP", "Basal",        1e8,  50.0, 250.0, 315.0),
    "rubisco_fix":    BiokineticReaction("Fijación CO2",  "Fotosíntesis", 1e10, 65.0, 200.0, 310.0),
    "nitrogenase":    BiokineticReaction("Fijación N2",   "Ciclo N",      1e12, 75.0, 300.0, 320.0),
    "metanogenesis":  BiokineticReaction("Metanogénesis", "Metabolismo anaerobio", 5e7, 45.0, 210.0, 305.0),
    "redox_sulfatos": BiokineticReaction("Sulfato Reductasa", "Respiración anaerobia", 2e8, 55.0, 230.0, 312.0),
    "fosforilacion_oxid": BiokineticReaction("ATP Sintasa", "Metabolismo aerobio", 1e9, 52.0, 260.0, 318.0),
    "ureasa": BiokineticReaction("Ureasa", "Catálisis extrema", 5e9, 48.0, 290.0, 335.0),
    "glucolitic ferment": BiokineticReaction("Hexocinasa", "Fermentación", 1e7, 42.0, 200.0, 300.0)
}



# 3. FUNCIONES BASE (RADIACIÓN, PROBABILIDAD Y CINÉTICA)

def planck_spectral_radiance(lam_meters: np.ndarray, T_star: float) -> np.ndarray:
    exponent = (H * C) / (lam_meters * K_B * T_star)
    return (2.0 * H * C**2 / (lam_meters**5)) / (np.exp(exponent) - 1.0)

def photon_energy_ev(lam_meters: np.ndarray) -> np.ndarray:
    return (H * C / lam_meters) / EV_TO_J

def single_bond_probability(E_ev: np.ndarray, E_b_ev: float, E_exc_ev: float = 1.8, 
                            beta1: float = 10.0, beta2: float = 10.0) -> np.ndarray:
    p_activation = 1.0 / (1.0 + np.exp(-beta1 * (E_ev - E_exc_ev)))
    p_survival = 1.0 / (1.0 + np.exp(beta2 * (E_ev - E_b_ev)))
    return p_activation * p_survival

def total_photochemical_window(E_ev: np.ndarray, bonds: Dict[str, ChemicalBond], 
                                E_exc_ev: float = 1.8) -> np.ndarray:
    prod_term = np.ones_like(E_ev, dtype=float)
    for bond in bonds.values():
        w_i = single_bond_probability(E_ev, bond.energy_ev, E_exc_ev)
        prod_term *= (1.0 - w_i)
    return 1.0 - prod_term

def eval_johnson_lewin(T_k: float, reaction: BiokineticReaction) -> float:
    E_a_j = reaction.E_a_kj_mol * 1000.0
    dH_denat_j = reaction.dH_denat_kj_mol * 1000.0
    C_factor = np.exp(dH_denat_j / (R_GAS * reaction.T_m_k))
    activation = reaction.A_pre_exp * np.exp(-E_a_j / (R_GAS * T_k))
    denaturation = 1.0 + C_factor * np.exp(-dH_denat_j / (R_GAS * T_k))
    return activation / denaturation



# 4. MÓDULOS DE SIMULACIÓN (HABITABILIDAD Y ÓRBITA)

class AstrobiologicalHabitabilityModel:
    def __init__(self, T_star: float, R_star: float, distance_au: float, 
                    R_planet: float, T_planet: float, albedo: float = 0.3):
        self.T_star = T_star
        self.R_star = R_star
        self.distance_m = distance_au * AU_M
        self.R_planet = R_planet
        self.T_planet = T_planet
        self.albedo = albedo
        self.A_cross_section = np.pi * (R_planet ** 2)
        self.A_surface = 4.0 * np.pi * (R_planet ** 2)

    def evaluate(self, bonds: Dict[str, ChemicalBond], reactions: Dict[str, BiokineticReaction] = {}, 
                    E_exc_ev: float = 1.8) -> dict:
        
        lam_grid = np.linspace(100e-9, 2000e-9, 1000)
        E_grid_ev = photon_energy_ev(lam_grid)
        b_lam = planck_spectral_radiance(lam_grid, self.T_star)
        
        flux_star_surface = np.pi * np.trapezoid(b_lam, lam_grid)
        flux_incidente = flux_star_surface * (self.R_star / self.distance_m)**2

        power_intercepted = flux_incidente * self.A_cross_section
        power_absorbed = power_intercepted * (1.0 - self.albedo)
        flux_surface_avg = power_absorbed / self.A_surface

        w_t_lam = total_photochemical_window(E_grid_ev, bonds, E_exc_ev)
        total_spectral_energy = np.trapezoid(b_lam, lam_grid)
        useful_spectral_energy = np.trapezoid(w_t_lam * b_lam, lam_grid)
        photochemical_index = float(useful_spectral_energy / total_spectral_energy)

        negentropy_rate_global = (4.0 / 3.0) * power_absorbed * ((1.0 / self.T_planet) - (1.0 / self.T_star))

        kinetic_rates = {key: eval_johnson_lewin(self.T_planet, rx) for key, rx in reactions.items()}

        return {
            "flux_incidente_W_m2": flux_incidente,
            "power_intercepted_GW": power_intercepted / 1e9,
            "power_absorbed_GW": power_absorbed / 1e9,
            "flux_surface_avg_W_m2": flux_surface_avg,
            "photochemical_index": photochemical_index,
            "negentropy_rate_global_MW_K": negentropy_rate_global / 1e6,
            "kinetic_rates": kinetic_rates
        }


class OrbitalKinematics:
    def __init__(self, semi_major_axis_m: float, eccentricity: float, 
                    rotation_period_s: float, orbital_period_s: float,
                    theta_rad: float, phi_rad: float, R_planet_m: float):
        self.a = semi_major_axis_m
        self.e = eccentricity
        self.T_rot = rotation_period_s
        self.T_orb = orbital_period_s
        self.theta = theta_rad
        self.phi = phi_rad
        self.R_p = R_planet_m
        
        v_raw = np.array([1.0 / np.tan(self.phi), 1.0, np.tan(self.theta)])
        self.v_hat = v_raw / np.linalg.norm(v_raw)

    def kepler_equation(self, E: float, M: float) -> float:
        return E - self.e * np.sin(E) - M

    def true_anomaly(self, t_orbit: float) -> float:
        M = (2 * np.pi / self.T_orb) * t_orbit
        E_initial = M if self.e < 0.8 else np.pi 
        E = fsolve(self.kepler_equation, E_initial, args=(M,))[0]
        tan_f_half = np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E / 2.0)
        return 2.0 * np.arctan(tan_f_half)

    def r_centro(self, t_orbit: float) -> np.ndarray:
        f = self.true_anomaly(t_orbit)
        r_t = (self.a * (1 - self.e**2)) / (1 + self.e * np.cos(f))
        return np.array([r_t * np.cos(f), r_t * np.sin(f), 0.0])

    def maximum_shadow_angle(self, num_samples: int = 1000) -> float:
        t_array = np.linspace(0, self.T_orb, num_samples)
        max_deviation = 0.0
        
        for t in t_array:
            R_c = self.r_centro(t)
            cos_psi = np.clip(np.dot(self.v_hat, -R_c / np.linalg.norm(R_c)), -1.0, 1.0)
            deviation = np.abs((np.pi / 2.0) - np.arccos(cos_psi))
            if deviation > max_deviation:
                max_deviation = deviation
                
        return max_deviation

    def permanent_darkness_area(self) -> dict:
        max_dev = self.maximum_shadow_angle()
        rho_perm = (np.pi / 2.0) - max_dev
        
        area = 0.0 if rho_perm < 0 else 2 * np.pi * (self.R_p**2) * (1.0 - np.cos(rho_perm))
        total_area = 4 * np.pi * (self.R_p**2)
        
        return {
            "darkness_area_m2": area,
            "percentage_darkness": (area / total_area) * 100,
        }



# 5. SISTEMA CENTRALIZADO Y EJECUCIÓN

def simular_planeta(config: PlanetConfig):
    print(f"\n==================================================")
    print(f"      EVALUACIÓN: {config.name.upper()} (T_med = {config.T_planet} K)")
    print(f"==================================================")
    
    # 1. Configurar Modelos
    model = AstrobiologicalHabitabilityModel(
        T_star=config.T_star, R_star=config.R_star_sun * R_SUN, distance_au=config.distance_au, 
        R_planet=config.R_planet_earth * R_EARTH, T_planet=config.T_planet, albedo=config.albedo
    )
    
    kin = OrbitalKinematics(
        semi_major_axis_m=config.distance_au * AU_M, eccentricity=config.eccentricity,
        rotation_period_s=config.rotation_period_s, orbital_period_s=config.orbital_period_s,
        theta_rad=np.deg2rad(config.theta_deg), phi_rad=np.deg2rad(config.phi_deg), 
        R_planet_m=config.R_planet_earth * R_EARTH
    )

    # 2. Ejecutar Cálculos (química CHONPS, con cinética biológica)
    resultados_hab = model.evaluate(BONDS_CHONPS, BIO_REACTIONS)
    # 2b. Evaluación paralela con química hipotética basada en silicio
    #     (mismo modelo estelar/planetario, solo cambia la tabla de enlaces)
    resultados_hab_si = model.evaluate(BONDS_SILICON)
    resultados_orb = kin.permanent_darkness_area()

    # 3. Presentar Datos
    print(f"Índice W_T Ponderado (CHONPS):   {resultados_hab['photochemical_index']:.4f}")
    print(f"Índice W_T Ponderado (Silicio):  {resultados_hab_si['photochemical_index']:.4f}")
    print(f"Potencia Absorbida:           {resultados_hab['power_absorbed_GW']:,.2f} GW")
    print(f"Generación Negentropía:       {resultados_hab['negentropy_rate_global_MW_K']:,.2f} MW/K")
    print(f"Área Oscuridad Permanente:    {resultados_orb['percentage_darkness']:.2f} %")
    print("-" * 50)
    
    for rx_name, rate in sorted(resultados_hab['kinetic_rates'].items(), key=lambda item: item[1]):
        print(f"Tasa r(T) [{rx_name}]: {rate:.2e} act/s")


if __name__ == "__main__":
    # Catálogo extensible de planetas
    planetas = [    # consultar los @dataclass si se tienen dudas sobre que indica cada param
        PlanetConfig("Tierra",      5778.0, 1.0,    1.0,    1.0,    288.0, 0.3,  0.0167, 86400,    3.154e7, 23.44, 90.0),
        PlanetConfig("Marte",       5778.0, 1.0,    1.52,   0.532,  214.0, 0.25, 0.0934, 88775,    5.935e7, 25.19, 90.0),
        PlanetConfig("Urano",       5778.0, 1.0,    19.23,  3.980,  68.00, 0.3, 0.044405, 62040,  2661041808, 97, 0),
        PlanetConfig("Kepler-452b", 5755.0, 1.1,    1.04,   1.63,   265.0, 0.3,  0.04,   86400*20, 3.32e7,  23.5,  90.0)
      # PlanetConfig("nombre", T_estrella, R_estrella_sol, dist_au, R_planet_tierra, T_planeta, albedo BOND, eccentr., periodo_rot_s, periodo_orbt_s, theta_grados, phi_grados.)
    ]

    for p in planetas:
        simular_planeta(p)