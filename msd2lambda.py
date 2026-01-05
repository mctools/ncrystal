import NCrystal as NC
import NCrystal.constants as nc_constants

def print_lambdas( cfgstr, temperature_k ):
    print(f'Investigating "{cfgstr}" at T={temperature_k}"')
    info = NC.createInfo(cfgstr+f';temp={temperature_k}K')
    k_B = nc_constants.constant_boltzmann
    hbar = nc_constants.constant_planck / nc_constants.k2Pi
    amu2eV = nc_constants.constant_dalton2eVc2 / nc_constants.constant_c**2
    for ai in info.atominfos:
        mass = ai.atomData.averageMassAMU()
        lambdaval = 2 * ai.msd * mass * amu2eV * k_B * temperature_k / hbar**2
        print(f'   lambda({ai.atomData.displayLabel()}) = {lambdaval:.4g}')


def main():
    T = 293.6
    cfgstrs = [
        'stdlib::C_sg194_pyrolytic_graphite.ncmat',
        'stdlib::Al_sg225.ncmat',
        'stdlib::Ag_sg225.ncmat',
        'stdlib::Pb_sg225.ncmat',
        'stdlib::Pd_sg225.ncmat',
        'stdlib::Sr_sg225.ncmat',
        'stdlib::Cu_sg225.ncmat',
        'stdlib::UO2_sg225_UraniumDioxide.ncmat',
    ]

    if True:
        #Extra UO2 data with VDOS tuned at 300K rather than 600K
        #(requires a plugin: "pip install ncrystal-plugin-UraniumOxideData")
        cfgstrs.append(
            'plugins::UraniumOxideData/UO2_sg225_UraniumOxide_vdos300K.ncmat'
        )

    for c in cfgstrs:
        print_lambdas( c, T )

if __name__=='__main__':
    main()
