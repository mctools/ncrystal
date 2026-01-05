import pandas as pd
import NCrystal as NC
from NCrystal.vdos import debyeIsotropicMSD
from scipy.constants import physical_constants, hbar, k

def get_constants():
    """
    Retrieves necessary physical constants using scipy.
    Returns a dictionary of constants in standard SI units,
    plus conversion factors for Angstroms and eV.
    """

    # 1 Dalton (AMU) in kg
    amu_kg = physical_constants['atomic mass constant'][0]

    # Electron volt in Joules (for energy conversion if needed,
    # though we stick to SI/Joule for the calculation to be safe)
    ev_joule = physical_constants['electron volt'][0]

    return {
        "hbar": hbar,          # J*s
        "k_B": k,              # J/K
        "amu_kg": amu_kg,      # kg
        "ev_joule": ev_joule   # J
    }

def msd_to_lambda(msd_angstrom2, mass_amu, temperature_k):
    """
    Converts Mean Square Displacement (MSD) <u^2> to the dimensionless
    phonon expansion coefficient lambda.

    Formula: lambda = MSD * (2 * M * k_B * T) / hbar^2

    Parameters:
    -----------
    msd_angstrom2 : float
        Mean Square Displacement in Angstroms squared (A^2).
    mass_amu : float
        Atomic mass in Atomic Mass Units (Daltons).
    temperature_k : float
        Temperature in Kelvin.

    Returns:
    --------
    lambda_val : float
        The dimensionless phonon expansion coefficient.
    """

    consts = get_constants()

    # Convert inputs to SI units (meters and kg)
    # 1 Angstrom = 1e-10 meters
    msd_m2 = msd_angstrom2 * (1e-10)**2
    mass_kg = mass_amu * consts["amu_kg"]

    # Calculate numerator: 2 * M * kB * T
    numerator = 2 * mass_kg * consts["k_B"] * temperature_k

    # Calculate denominator: hbar^2
    denominator = consts["hbar"]**2

    # Calculate lambda (units of m^2 * kg * J/K * K / (J*s)^2 ) -> dimensionless
    lambda_val = msd_m2 * (numerator / denominator)

    return lambda_val

def get_lambda_crystal(mat, atom_lbl, T):
    # Get atom information:
    atom_info = mat.info.findAtomInfo(atom_lbl)
    mass = atom_info.atomData.averageMassAMU()

    # Get the msd form
    #ORIG: debyeTemp = mat.info.atominfos[0].debyeTemperature
    debyeTemp = atom_info.debyeTemperature
    iso_msd = debyeIsotropicMSD(debye_temperature=debyeTemp, temperature=T, mass=mass)

    # Get the lambda:
    return msd_to_lambda(iso_msd, mass, T)

def get_lambda(base_ncmat_file, atom_lbl, T):
    lambda_test = {}
    for i in range(len(base_ncmat_file)):
        ncmat_filename = "".join([f'{base_ncmat_file[i]};', f'temp={T:.1f}K'])
        # NC: Loaded Material:
        mat = NC.load(ncmat_filename)

        # Get lambda from Ncrystal:
        # Note: 'solid_cinel' function was missing in snippet, assumed external
        key_name = "".join([atom_lbl[i], ''])
        lambda_test[key_name] = {
            'ncrystal': get_lambda_crystal(mat, atom_lbl[i], T),
            # 'solid_cinel': get_lambda_solid_cinel(mat, atom_lbl[i], T)
        }

    return pd.DataFrame(lambda_test)

base_ncmat_files = [
    'C_sg194_pyrolytic_graphite.ncmat',
    'Al_sg225.ncmat',
    'Ag_sg225.ncmat',
    'Pb_sg225.ncmat',
    'Pd_sg225.ncmat',
    'Sr_sg225.ncmat',
    'Cu_sg225.ncmat',
    'UO2_sg225_UraniumDioxide.ncmat',
    'UO2_sg225_UraniumDioxide.ncmat'
]

T = 293.6 # K
atom_lbls = ['C','Al', 'Ag', 'Pb', 'Pd', 'Sr', 'Cu', 'O', 'U']

print( get_lambda(base_ncmat_files,atom_lbls,T) )
