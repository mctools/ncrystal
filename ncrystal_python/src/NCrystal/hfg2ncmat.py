
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

"""

Utilities for generating NCMAT files for hydrogen rich amorphous materials
(based on doi:10.1088/1361-648X/abfc13).

"""

__all__ = [ 'hfg2ncmat']

def _default_debye_temp():
    #Note: keep in synch with value mentioned in docstring of hfg2ncmat
    #function, as well as the NCMATComposer.from_hfg method.
    return 400.0

def _parseChemicalFormula(chemform):
    import re
    def err():
        from .exceptions import NCBadInput
        raise NCBadInput(f'Invalid chemical formula: "{chemform}"')
    def is_ascii(s):
        #str.isascii requires python 3.7
        return all(ord(c) < 128 for c in s)
    if not is_ascii(chemform) or not chemform.isalnum():
        err()
    d={}
    prev=None
    def add(element,count):
        d[element] = d.get(element,0) + count

    chemform_unused = chemform.strip()
    for e in re.findall('[A-Z][a-z]?|[0-9]+', chemform.strip()):
        if chemform_unused.startswith(e):
            chemform_unused = chemform_unused[len(e):].strip()
        else:
            err()
        if e.isdigit():
            if prev is None or prev.isdigit():
                err()
            v=int(e)
            if v<1 or v>1000000000:
                err()
            add(prev,v)
            prev=None
            continue
        if prev is not None:
            add(prev,1)
            prev = None
        prev = e
    if prev is not None:
        add(prev,1)
    if chemform_unused:
        err()
    return d

def _parseHydrogenBindingSpecification(spec):
    from .exceptions import NCBadInput
    def err():
        raise NCBadInput(f'Invalid hydrogen binding specification: "{spec}"')
    def is_ascii(s):
        #str.isascii requires python 3.7
        return all(ord(c) < 128 for c in s)
    if not is_ascii(spec) or not spec.replace('+','').isalnum():
        err()

    fg2nH=dict(CHali=1,CHaro=1,CH2=2,CH3=3,NH=1,NH2=2,NH3=3,OH=1,SH=1)
    d={}
    for part in spec.split('+'):
        p=part.strip().split('x')
        if len(p)!=2:
            err()
        count,fg=p
        if fg not in fg2nH:
            raise NCBadInput(f'Unknown functional group in "{spec}": "{fg}"')
        if not count.isdigit() or not 1 <= int(count) <= 10000:
            raise NCBadInput(f'Invalid count in "{spec}": "{count}"')
        d[fg] = d.get(fg,0)+int(count)*fg2nH[fg]
    return d

def hfg2ncmat( spec,
               formula, *,
               density,
               title,
               debyetemp = _default_debye_temp(),
               verbose = True,
               notrim = False ):
    """Function which can be used to generate NCMAT data for hydrogen-rich
amorphous materials, in which the hydrogen atoms are bound to certain standard
functional groups (e.g. carbohydrates, polyimides, polymers, ...). Based on the
material's density, (empirical) chemical formula, and the specification of
hydrogen bindings in terms of standard functional groups, NCMAT data is
generated. In this NCMAT data, non-hydrogen atoms are treated with a simplistic
model (idealised Debye model of phonon vibrations, assuming a Debye temperature
of 400K for all atoms unless the debyetemp parameter is provided. The hydrogen
atoms are treated with a proper phonon density of state (VDOS) curve, which is
constructed based on the provided binding specifications. This is done using an
idea (and VDOS curves from) the following publication:

  "Thermal neutron cross sections of amino acids from average contributions
  of functional groups", G. Romanelli, et. al., J. Phys.: Condens. Matter,
  (2021). doi:10.1088/1361-648X/abfc13

Example of valid spec strings for the spec parameter. Note only the hydrogen
bindings are important here:

     "1xCHali+2xCH2+8xCHaro+2xCH3+1xOH" (20 hydrogen atoms)
     "1xNH+3xSH+2xCH2" (8 hydrogen atoms).

List of valid bindings which can be used in the spec string:

      CHali (1 hydrogen atom, aliphatic binding)
      CHaro (1 hydrogen atom, aromatic binding, e.g. on phenyl group).
      CH2 (2 hydrogen atoms)
      CH3 (3 hydrogen atoms)
      NH (1 hydrogen atom)
      NH2 (2 hydrogen atoms)
      NH3 (3 hydrogen atoms)
      OH (1 hydrogen atom)
      SH (1 hydrogen atom)

Note that as only hydrogen atoms will have a realistic VDOS curve, and as the
incoherent approximation is employed, the realism of the resulting material
modelling will be higher for materials with more hydrogen atoms. As a metric of
this, the function prints out how much of the total scattering cross section
(sum of sigma_coh+sigma_incoh for all atoms), is due to the sigma_incoh
contribution from hydrogen atoms (unless verbose=False).

"""
    #NOTE: Keep the description above synchronised with help text of of the
    #ncrystal_hfg2ncmat command-line app!

    from ._hfgdata import get_data as hfg_get_data

    if verbose:
        from ._common import print
    else:
        def print(*args,**kwargs):
            pass

    spec_parsed = _parseHydrogenBindingSpecification(spec)
    chemform = _parseChemicalFormula(formula)

    nH = chemform.get('H',-999)
    if sum(v for k,v in spec_parsed.items()) != nH:
        from .exceptions import NCBadInput
        raise NCBadInput('number of hydrogen items in --formula'
                         ' and --spec does not match up.')

    egrid,vdos = hfg_get_data()

    #Calculate hydrogen VDOS:

    h_vdos = sum( nh * vdos[fgname] for fgname,nh in sorted(spec_parsed.items()))
    h_vdos /= h_vdos.max()

    debye_cut_bin = 4#shave off first few points, to force parabolic behaviour near E=0.
    if not notrim:
        for i in range(len(h_vdos)):
            if egrid[i]<0.05:
                continue#don't touch first bins
            if h_vdos[i]<1e-7:
                h_vdos[i]=0.0
        ntrim=2
        assert debye_cut_bin%ntrim==0
        debye_cut_bin //= ntrim
        h_vdos = h_vdos[::ntrim]
        egrid = egrid[::ntrim]

    while len(h_vdos)>10 and h_vdos[-2]==0.0 and h_vdos[-1]==0.0:
        h_vdos = h_vdos[0:-1]
        egrid = egrid[0:-1]

    appname='ncrystal_hfg2ncmat'

    initial_comments='#\n'
    for ll in title.replace('\\n','\n').splitlines():
        initial_comments+= f'# {ll}\n' if ll.strip() else '#\n'
    initial_comments+='#'

    recordargs=[ f'--spec={spec}',
                 f'--formula={formula}',
                 f'--density={density}' ]
    if _default_debye_temp() != debyetemp:
        recordargs.append('--debyetemp=%.14g'%debyetemp)
    if notrim:
        recordargs.append('--notrim')
    recordargs='\n#  '.join(recordargs)

    out = f"""NCMAT v5
{initial_comments}
# ----------------------------------------------------------------------------
#
# File generated by {appname} with:
#
#  {recordargs}
#
# The hydrogen VDOS curve in this file is constructed using data from:
#
#   "Thermal neutron cross sections of amino acids from average contributions
#   of functional groups", G. Romanelli, et. al., J. Phys.: Condens. Matter,
#   (2021). doi:10.1088/1361-648X/abfc13
#
# Considering only scat. lengths and composition, incoherent scattering on
# hydrogen accounts for <<<HYDROGENPERCENTAGE>>> of scattering in this material (a higher value
# implies more accurate modelling).
#
@STATEOFMATTER
  solid
@DENSITY
  {density} g_per_cm3
"""

    ntot = sum(n for elem,n in chemform.items())
    for elem,n in chemform.items():
        if elem=='H':
            continue
        out+=f"""@DYNINFO
  element {elem}
  fraction {n}/{ntot}
  type vdosdebye
  debye_temp {debyetemp}
"""
    out+=f"""@DYNINFO
  element  H
  fraction {nH}/{ntot}
  type     vdos
  vdos_egrid {0.001*egrid[debye_cut_bin]:.8g} {0.001*egrid[-1]:.8g}
"""
    from .ncmat import formatVectorForNCMAT
    out += formatVectorForNCMAT('vdos_density',h_vdos[debye_cut_bin:])
    print ("Generated NCMAT data...")
    print ("Verifying that it can be loaded with NCrystal...")
    from .core import load as nc_load
    loaded_mat = nc_load(out)
    print ("Succesfully loaded...")

    breakdown = []
    for di in loaded_mat.info.dyninfos:
        for cohtype in 'incoherent','coherent':
            key='%s (%s)'%(di.atomData.displayLabel(),cohtype)
            sigma = di.atomData.coherentXS() if cohtype=='coherent' else di.atomData.incoherentXS()
            breakdown+=[(di.fraction*sigma,key)]
    sigma_tot = sum(s for s,k in breakdown)
    print("Contribution breakdown based on composition:")
    for s,key in reversed(sorted(breakdown)):
        print ("  Contribution to bound scattering XS from %s is %5.2f %%"%(key.ljust(14),s*100.0/sigma_tot))
        if key=='H (incoherent)':
            contrib_incH = s/sigma_tot

    return out.replace('<<<HYDROGENPERCENTAGE>>>','%.3g%%'%(100.0*contrib_incH))
