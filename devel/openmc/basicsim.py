
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

import numpy as np
import openmc

TRUE_VALUES = np.array([5.60481868e+01, 5.75019668e-04, 7.67958251e-01])
TRUE_ERRORS = np.array([2.58956360e-03, 4.83597510e-06, 1.96171839e-03])
TRUE_POSITIONS = np.array([0.00872665, 1.15202318, 2.29531971])

def is_interactive():
    import sys
    return ( '--interactive' in sys.argv[1:]
             or '-i' in sys.argv[1:] )

def pencil_beam_model(cfg, E0, N):
    """Return an openmc.Model() object for a monoenergetic pencil
     beam hitting a 1 mm sphere filled with the material defined by
     the cfg string, and compute angular distribution for the 
     Debye-Scherrer cones"""

    import NCrystal as NC

    # Material definition

    m1 = openmc.Material.from_ncrystal(cfg)
    materials = openmc.Materials([m1])

    # Geometry definition

    sample_sphere = openmc.Sphere(r=0.1)
    outer_sphere = openmc.Sphere(r=100, boundary_type="vacuum")
    cell1 = openmc.Cell(region=-sample_sphere, fill=m1)
    cell2_region = +sample_sphere & -outer_sphere
    cell2 = openmc.Cell(region=cell2_region, fill=None)
    geometry = openmc.Geometry([cell1, cell2])

    # Source definition

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, -20))
    source.angle = openmc.stats.Monodirectional(reference_uvw=(0, 0, 1))
    source.energy = openmc.stats.Discrete([E0], [1.0])

    # Execution settings

    settings = openmc.Settings()
    settings.source = source
    settings.run_mode = "fixed source"
    settings.batches = 10
    settings.particles = N

    # Tally definition

    wl = NC.ekin2wl(E0)
    angular_binwdidth = np.radians(1.0)
    debye_cones = np.array([ 2.0*np.asin(wl/(2.0*hkl.d)) for hkl in NC.createInfo(cfg).hklObjects() if wl < 2.0*hkl.d ])
    angular_bins = np.sort(np.concatenate(([0,angular_binwdidth], debye_cones-0.5*angular_binwdidth, debye_cones+0.5*angular_binwdidth)))
    
    tally1 = openmc.Tally(name="debye-scherrer cones")
    tally1.scores = ["current"]
    filter1 = openmc.SurfaceFilter(sample_sphere)
    filter2 = openmc.PolarFilter(angular_bins)
    filter3 = openmc.CellFromFilter(cell1)
    tally1.filters = [filter1, filter2, filter3]

    tally2 = openmc.Tally(name="angular distribution")
    tally2.scores = ["current"]
    filter2b = openmc.PolarFilter(np.linspace(0, np.pi, 180+1))
    tally2.filters = [filter1, filter2b, filter3]
    tallies = openmc.Tallies([tally1, tally2])

    return openmc.Model(geometry, materials, settings, tallies)

def check_results(sp_file):
    with openmc.StatePoint(sp_file) as sp:
        tal = sp.get_tally(name='debye-scherrer cones')
        df = tal.get_pandas_dataframe()
    bin_high = df['polar high [rad]'].values
    bin_low = df['polar low [rad]'].values
    values = df['mean'].values/(bin_high - bin_low)
    errors = df['std. dev.'].values/(bin_high - bin_low)
    assert np.shape(TRUE_VALUES) == np.shape(values), 'Wrong number of Debye-Scherrer cones'
    test_result = np.all(values > TRUE_VALUES-TRUE_ERRORS) and np.all(values < TRUE_VALUES+TRUE_ERRORS)
    if test_result:
        print('OpenMC run succesful.')
    else:
        print('OpenMC run not succesful.')
    print('Computed intensity for Debye-Scherrer cones:')
    for true_pos,true_val,true_err,val,err in \
            zip(TRUE_POSITIONS, TRUE_VALUES, TRUE_ERRORS, values, errors):
        print(f' Cone at {np.rad2deg(true_pos):>5.1f} deg / Expected: {true_val:.2e} +/- {true_err:.2e} / Computed: : {val:.2e} +/- {err:.2e}')
    return test_result 

def plot_tally(sp_file):
    with openmc.StatePoint(sp_file) as sp:
        tal = sp.get_tally(name='debye-scherrer cones')
        df = tal.get_pandas_dataframe()
    bin_high = df['polar high [rad]'].values
    bin_low = df['polar low [rad]'].values
    cone_positions = 0.5*(bin_high + bin_low)
    cone_values = df['mean'].values/(bin_high - bin_low)
    cone_errors = df['std. dev.'].values/(bin_high - bin_low)

    with openmc.StatePoint(sp_file) as sp:
        tal = sp.get_tally(name='angular distribution')
        df = tal.get_pandas_dataframe()
    bin_high = df['polar high [rad]'].values
    bin_low = df['polar low [rad]'].values
    angdist_bins = bin_high
    angdist_values = df['mean'].values/(bin_high - bin_low)
    angdist_errors = df['std. dev.'].values/(bin_high - bin_low)
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.step(np.rad2deg(bin_high), angdist_values)
    plt.errorbar(np.rad2deg(cone_positions), cone_values, yerr=cone_errors, fmt='.')
    plt.errorbar(np.rad2deg(TRUE_POSITIONS), TRUE_VALUES, yerr=TRUE_ERRORS, fmt='.', label='True values')
    plt.yscale('log')
    plt.ylabel('density')
    plt.xlabel('Polar angle [deg]')
    plt.legend()
    plt.show()

def main():
    
    n_particles = 1000000
    E0 = 0.0045  # eV
    cfg = 'Al_sg225.ncmat'
    model = pencil_beam_model(cfg, E0, n_particles)
    sp_file = model.run()
    
    if is_interactive():
        plot_tally(sp_file)
        
    assert check_results(sp_file), 'Statistically significant deviation from expected values'

if __name__ == '__main__':
    main()
