
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

def main():

    import openmc

    mat_al = openmc.Material.from_ncrystal('Al_sg225.ncmat')
    materials = openmc.Materials([mat_al])

    sphere0 = openmc.Sphere(r=10, boundary_type='vacuum')
    cell0 = openmc.Cell(region=-sphere0, fill=mat_al)
    geometry = openmc.Geometry([cell0])

    settings = openmc.Settings()
    settings.source = openmc.IndependentSource()
    settings.run_mode = 'fixed source'
    settings.batches = 10
    settings.particles = 1000

    # Run OpenMC model
    model = openmc.Model(geometry, materials, settings)
    model.run()

if __name__ == '__main__':
    main()
