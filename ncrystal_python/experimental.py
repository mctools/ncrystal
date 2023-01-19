"""Experimental utility code"""


################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2022 NCrystal developers                                   ##
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


#Import NCrystal module with importlib to make sure we always get the one from
#the same directory as the current file:
import pathlib
import math
import copy
from collections import Counter as collections_Counter

#Main NCrystal module
import importlib
__name = 'NCrystal.experimental' if __name__=='__main__' else __name__
NC = importlib.import_module( '..', __name )#importing __init__.py
nc_common = importlib.import_module( '.._common', __name )
print = nc_common.print

__all_known_element_names=[None,None]#list,set
def knownElementNames():
    """Returns tuple of all known element names ordered by increasing Z value."""
    if __all_known_element_names[0] is None:
        l=[]
        for a in NC.iterateAtomDB():
            if a.isNaturalElement():
                l.append( (a.Z(),a.elementName()) )
        l = tuple( name for z,name in sorted(l) )
        __all_known_element_names[0] = l
        __all_known_element_names[1] = frozenset(l)
    return __all_known_element_names[0]

def isKnownElement( label ):
    if __all_known_element_names[1] is None:
        knownElementNames()#trigger init
    return label in __all_known_element_names[1]

__all_element_names=tuple(["H", "He", "Li", "Be", "B", "C",
                           "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
                           "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
                           "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
                           "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
                           "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                           "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
                           "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                           "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
                           "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
                           "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
                           "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts",
                           "Og"])
def allElementNames():
    """Like knownElementNames, but also returns names of elements for which NCrystal has no data."""
    return __all_element_names

class NCMATComposer:
    #TODO: embedded cfg?

    def __init__(self, params = None ):
        self._ichange = 0#increment on all changes, for downstream caching
        self.__params = copy.deepcopy(params) if params else {}
        self.__loadcache = None

    def clone(self):
        o = NCMATComposer( params = self.__params )
        o._ichange = self._ichange
        o.__loadcache = self.__loadcache
        return o

    def __dirty(self):
        self._ichange += 1
        self.__loadcache = None

    def get_cache_key( self ):
        """
        Returns tuple of integers which will be unchanged if generated NCMAT
        content is unchanged
        """
        return ( int(id(self), self._ichange) )

    def to_dict(self):
        return copy.deepcopy(self.__params)

    def set_cellsg( self, *, a,b,c, alpha,beta,gamma, spacegroup=None ):
        #TODO: input validation! Also ensure high precision in numbers below
        """Assumes alpha,beta,gamma in degrees, a,b,c in Aangstrom, spacegroup an int in range 1..230."""
        self.__dirty()
        self.__params['cellsg'] = dict(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,spacegroup=spacegroup)

    def set_cellsg_cubic( self, a, *, spacegroup=None ):
        assert spacegroup is None or 195<=int(spacegroup)<=230
        self.set_cellsg( a=a,b=a,c=a,alpha=90.,beta=90.,gamma=90.,
                         spacegroup = ( int(spacegroup) if spacegroup else None ) )

    def __lines_cellsg(self):
        c = self.__get_cellsg()
        l = f"""
        @CELL
        lengths {c['a']:g} {c['b']:g} {c['c']:g}
        angles {c['alpha']:g} {c['beta']:g} {c['gamma']:g}
        """
        if c['spacegroup']:
            l+=f"""
            @SPACEGROUP
            {c['spacegroup']}
            """
        return l

    def __add_dyninfo( self, atomlabel, dyninfo ):
        self.__dirty()
        if not 'dyninfos' in self.__params:
            self.__params['dyninfos'] = { atomlabel : dyninfo }
        else:
            self.__params['dyninfos'][atomlabel] = dyninfo

    #def set_dyninfo_sterile( self, atomlabel ):
    #    self.__dyninfos[atomlabel]=dict(ditype='sterile')

    #TODO: def set_dyninfo_scatknl( self, atomlabel ):
    #    self.__dyninfos[atomlabel]=dict(ditype='scatknl')

    def set_dyninfo_vdos( self, atomlabel, vdos_egrid, vdos, comment = None ):
        self.__add_dyninfo( atomlabel, dict( ditype='vdos',
                                             vdos_egrid=vdos_egrid,
                                             vdos=vdos,
                                             comment = None if not comment else str(comment) ) )

    def set_dyninfo_vdosdebye( self, atomlabel, debye_temp, comment = None ):
        self.__add_dyninfo( atomlabel, dict( ditype='vdosdebye',
                                             debye_temp = debye_temp,
                                             comment = None if not comment else str(comment) ) )

    def set_dyninfo_msd( self, atomlabel, *, msd, temperature, mass=None, comment = None ):
        """Calculate and set Debye temperature based on
        mean-squared-displacement (msd) value (in Aa^2). This also needs the temperature
        value for which the msd value is associated (in kelvin), as well as the atom mass
        (the mass parameter can be omitted in case the atom label is the name of
        an known element of isotope like H, D, O16, ...).
        """
        if mass is None:
            atomdata = NC.atomDB( atomlabel, throwOnErrors=False )
            if atomdata is None:
                raise ValueError( 'The set_dyninfo_msd method requires a mass (amu) value'
                                  + ' when the atomlabel is not a known element or isotope name' )
            mass = atomdata.averageMassAMU()

        dt = NC.debyeTempFromIsotropicMSD( msd = msd, mass = mass,
                                           temperature = temperature )

        c = f'Debye temperature value derived from msd={msd:g}Aa^2 @ T={temperature:g}K (mass={mass:g}u)'
        if comment:
            c=f'{comment.strip()}\n{c}'
        self.set_dyninfo_vdosdebye( atomlabel, dt, comment = c )

    def set_atompos( self, atompos ):

        """The atompos parameter must be list of entries
        (atomlabel,x,y,z) or (atomlabel,x,y,z,site_occupancy).

        If provided, the site_occupancy value (a number in (0,1] must be the
        same for all identical atomlabels.

        The site_occupancy parameter is for now highly experimental and will
        result in generated NCMAT files with fake sterile atoms inserted for
        technical reasons (resulting in a higher density but hopefully correct
        physics as long as the user does not subsequently try to override the
        density directly via the cfg-string parameters).
        """

        pos,occumap = [],{}
        def remap_pos(x):
            x = float(x)
            if -1.0<=x<0.0:
                x += 1.0
            if 1.0<=x<2.0:
                x -= 1.0
            return x

        for e in atompos:
            assert len(e) in (4,5)
            lbl,x,y,z = str(e[0]),remap_pos(e[1]),remap_pos(e[2]),remap_pos(e[3])

            if -1.0<=x<0.0:
                x+=1.0
            pos.append( (lbl,x,y,z) )
            occu = float(e[4]) if len(e)==5 else 1.0
            if not ( 0.0 < occu <= 1.0 ):
                raise ValueError('site_occupancy values must be in (0.0,1.0]')
            if not lbl in occumap:
                occumap[lbl]=occu
            else:
                if not occumap[lbl] == occu:
                    raise ValueError('All site occupancies for atomlabel "%s" are not identical'%lbl)

        #remove occu==1.0 entries from map (and sort):
        occumap = dict( sorted( (k,v) for k,v in occumap.items() if v!=1.0 ) )
        if occumap:
            nc_common.warn('Support for site_occupancy is highly experimental. Do *not* attempt to directly override'
                           +' value with the "density" cfg-parameter for this material (unless it is to simply scale it)')
        self.__dirty()
        self.__params['atompos'] = dict( pos=pos, occumap=occumap )

    def __get_atompos(self):
        atompos = self.__params.get('atompos',None)
        if not atompos:
            raise NC.NCMissingInfo("Atom positions in cell missing (add via call to set_atompos)")
        return atompos

    #    def set_temperature(self,temperature): #should be for @TEMPERATURE section!
    #        """set temperature value in kelvins"""
    #        self.append_cfg(f';temp={temperature}')

#    def append_cfg(self,cfg):
#        self.__dirty()
#        existing = self.__params.get('cfgparams','')
#        if not existing:
#            self.__params['cfgparams'] = cfg
#        else
#            self.__params['cfgparams'] = '%s;%s'%(existing,cfg)
#
#    def __get_cfgparams(self):
#        return self.__params.get('cfgparams','')

    def __get_cellsg(self):
        cellsg = self.__params.get('cellsg',None)
        if not cellsg:
            raise NC.NCMissingInfo("Cell/spacegroup information missing (add via call to set_cellsg)")
        return cellsg

    def __get_dyninfos(self):
        dyninfos = self.__params.get('dyninfos',None)
        if not dyninfos:
            raise NC.NCMissingInfo("Dynamics information missing (add via calls to set_dyninfo_...")
        return dyninfos

    def __determine_labels_and_atomdb( self ):
        atompos = self.__get_atompos()
        used_labels = set( lbl for lbl,x,y,z in atompos['pos'] )
        occumap = atompos['occumap']
        lbl_map = {}
        #Anything with a site_occupancy!=1.0 needs a special label (X1,X2,...)
        #so we can emulate the site occupancies via sterile atoms in the @ATOMDB section.
        for lbl,occu in occumap.items():
            assert 0.0 < occu < 1.0
            for j in range(1,100+1):
                if j==100:
                    raise RuntimeError('Ran out of custom labels (X1..X99 all used).')
                x = 'X%i'%j
                if not x in used_labels:
                    lbl_map[lbl] = x
                    used_labels.add(x)
                    break
        #Find elements to use for the role of "sterile atoms":
        sterile_atom_map = {}
        for lbl,occu in occumap.items():
            if not isKnownElement( lbl ):
                #must be pure element for now, since we need to be able to look up the mass.
                raise RuntimeError('site_occupancies != 1.0 are currently implemented via a hack that only works for atom '
                                   +'labels that are pure known natural elements (offending label: "%s").'%lbl)
            #Hijack the most obscure (i.e. highest Z) unused element name for this hack:
            hijacked_elem = None
            for elem in allElementNames()[::-1]:
                if not elem in used_labels:
                    hijacked_elem = elem
                    used_labels.add( elem )
                    break
            if not hijacked_elem:
                raise RuntimeError('Could not find unused element name with which to implement hack for site_occupancy support')
            sterile_atom_map[ lbl] = hijacked_elem
        atomdb_lines = []
        for elem,hijacked_elem in sorted(sterile_atom_map.items()):
            atomdata = NC.atomDB(elem)
            mass = atomdata.averageMassAMU()
            occu = occumap[elem]
            atomdb_lines.append( f'#Emulating site_occupancy={occu:g} for {elem} by hijacking {hijacked_elem} to play the role as "sterile {elem}":' )
            atomdb_lines.append( f'{hijacked_elem} {mass:g}u 0fm 0b 0b # sterile but same mass as {elem}' )
            atomdb_lines.append( f'{lbl_map[elem]} is {occu:g} {elem} {1.0-occu:g} {hijacked_elem}' )
        return lbl_map,atomdb_lines

    def __lines_atompos(self,lbl_map):
        atompos = self.__get_atompos()
        l="@ATOMPOSITIONS\n"
        fmt = nc_common.prettyFmtValue
        for lbl,x,y,z in atompos['pos']:
            l += f'{lbl_map.get(lbl,lbl)} {fmt(x)} {fmt(y)} {fmt(z)}\n'
        return l

    def __lines_dyninfo(self,lbl_map):
        #fixme: check consistent atompos and dyninfo labels (i.e. no unused dyninfos added).
        #base fractions on atompos entries:
        atompos = self.__get_atompos()
        dyninfos = self.__get_dyninfos()
        natoms = len(atompos['pos'])
        counts = collections_Counter(list(e[0] for e in atompos['pos']))
        lines=''
        for lbl,count in sorted(counts.items()):
            frac = float(count) / natoms
            lines += '@DYNINFO\n'
            dyninfo = dyninfos.get(lbl,None)
            if not dyninfo:
                raise NC.NCMissingInfo(f'Missing dyninfo for atom label "{lbl}" (add via call to set_dyninfo_...)')
            comment = dyninfo['comment']
            if comment:
                for e in comment.splitlines():
                    lines += f'  # {e.strip()}\n'

            lines += f'element {lbl_map.get(lbl,lbl)}\n'
            gcd = math.gcd(count,natoms)
            if gcd == natoms:
                lines += 'fraction 1\n'
            else:
                lines += 'fraction %i/%i\n'%(count//gcd,natoms//gcd)
            ditype = dyninfo['ditype']
            lines += f'type {ditype}\n'
            if ditype == 'vdosdebye':
                dt=dyninfo['debye_temp']
                lines += f'debye_temp {dt:g}\n'
            elif ditype == 'vdos':
                assert False,"not implemented"
                pass#TODO
#  vdos_egrid   0.013300161256 0.38127128934268
#  vdos_density 0.2594369576 0.3556097936 0.4375132376 0.4857022268 0.4915645468
#               0.444171744 0.340902404 0.2143757752 0.1140493096 0.0607265688
#               <...>
#               0.0593700356 0.0693096064 0.0812168444 0.0885391804 0.082003864
#               0.0920026688 0.0224714976 0.0040273484
            elif ditype == 'scatknl':
                assert False,"not implemented"
        return lines

    def register_as( self, virtual_filename ):
        NC.registerInMemoryFileData( virtual_filename, self.create_ncmat() )

    def load(self,cfg_param='',force = False):
        #cfg_param = ";dcutoff=0.5;comp=elas
        if force or self.__loadcache is None or self.__loadcache[0] != cfg_param:
            self.__loadcache = ( cfg_param, NC.directMultiCreate( self.create_ncmat(), cfg_param ) )
        return self.__loadcache[1]

    def create_ncmat(self):
        lbl_map,atomdb_lines = self.__determine_labels_and_atomdb()

        l  = self.__lines_cellsg()
        l += self.__lines_atompos( lbl_map )
        l += self.__lines_dyninfo( lbl_map )
        if atomdb_lines:
            l += '@ATOMDB\n'
            l += '\n'.join(atomdb_lines)
            l += '\n'

        out=["NCMAT v5\n#Autogenerated by %s"%self.__class__.__name__]
        for e in l.splitlines():
            e=e.strip()
            if e.startswith('@'):
                out.append(e)
            elif e:
                out.append('  '+e)
        return '\n'.join(out)+'\n'

#    def load(self,cfg_param=''):
#        #cfg_param = ";dcutoff=0.5;comp=elas
#        if self.__ncobj is not None and self.__ncobj[0] == cfg_param:
#            return self.__ncobj[1]
#            self.__ncobj = ( cfg_param, NC.directMultiCreate( self.create_ncmat(), cfg_param ) )
#        return self.__ncobj

class NCHelper_SimpleCrystalMaterial:
    def __init__(self):
        self.__ncmat = NCMATComposer()
        self.__temperature = 293.15
        self.__dcutoff = None
        self.__usrcfgparam = None
        #enable.disable absorption

    def clone(self):
        o = NCHelper_SimpleCrystalMaterial()
        o.__ncmat = self.__ncmat.clone()
        o.__temperature = self.__temperature
        o.__dcutoff = self.__dcutoff
        o.__usrcfgparam = self.__usrcfgparam
        return o

    def set_cellsg( self, *, a,b,c, alpha,beta,gamma, spacegroup=None ):
        """Assumes alpha,beta,gamma in degrees, a,b,c in Aangstrom, spacegroup an int in range 1..230."""
        self.__ncmat.set_cellsg(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,spacegroup=spacegroup)

    def set_cellsg_cubic( self, a, *, spacegroup=None ):
        self.__ncmat.set_cellsg_cubic( a=a, spacegroup=spacegroup )

    def set_dyninfo_vdos( self, atomlabel, vdos_egrid, vdos ):
        self.__ncmat.set_dyninfo_vdos( atomlabel, vdos_egrid, vdos )

    def set_dyninfo_vdosdebye( self, atomlabel, debye_temp ):
        self.__ncmat.set_dyninfo_vdosdebye( atomlabel, debye_temp )

    def set_atompos( self, atompos ):
        self.__ncmat.set_atompos( atompos )

    def set_temperature( self, temperature ):
        self.__temperature = temperature

    def set_dcutoff( self, dcutoff ):
        self.__dcutoff = dcutoff

    def set_mean_squared_displacement( self, atomlabel, msd, mass=None, temperature=None ):
        if not isKnownElement( atomlabel ):
            #NB: We could perhaps eventually relax this:
            raise ValueError('The set_mean_squared_displacement function requires a mass (amu) value when the atomlabel is not a known element name')

        if mass is None:
            atomdata = NC.atomDB(atomlabel)
            mass = atomdata.averageMassAMU()
        dt = NC.debyeTempFromIsotropicMSD( msd = msd,
                                           temperature = self.__temperature if temperature is None else temperature,
                                           mass = mass )
        self.set_dyninfo_vdosdebye( atomlabel, dt )

    def _cfgstr(self):
        l=[]
        if self.__usrcfgparam is not None:
            l.append(self.__usrcfgparam)
        if self.__dcutoff is not None:
            l.append(f'dcutoff={self.__dcutoff}')
        l.append(f'temp={self.__temperature}')
        return ';'.join(l)

    def load( self ):
        return self.__ncmat.load(self._cfgstr())

    @property
    def info( self ):
        return self.load().info

    def scatter(self, *, ekin=None,direction=None,wl=None,repeat=None):
        return load().scatter.scatter(ekin=ekin,direction=direction,wl=wl,repeat=repeat)

    def xsect( self, *, ekin=None, direction=None, wl=None, repeat=None, macro=False, absn=True, scat=True ):
        #Macroscopic as well!
        o = self.load()
        kwargs=dict(ekin=ekin, direction=direction, wl=wl, repeat=repeat)
        xsa = o.absorption.xsect(**kwargs) if absn else None
        xss = o.scatter.xsect(**kwargs) if scat else None
        if xsa is None:
            xs = xss
        elif xss is None:
            xs = xsa
        else:
            xs = xss + xsa
        if xs is None:
            return
        if not macro:
            return xs
        nd = o.info.getNumberDensity()
        return xs / nd if nd > 0.0 else ( xs + math.infinity )

def _main():
    #Fixme: to unit test
    x = NCMATComposer()
    x.set_cellsg_cubic(4.0,spacegroup=225)
    x.set_atompos( [ ('Al',   0, 0.5, 0.5),
                     ('Al',   0,   0,   0),
                     ('Al', 0.5, 0.5,   0),
                     ('Al', 0.5,   0, 0.5) ] )
    #x.set_dyninfo_freegas('Al')
    #x.set_dyninfo_vdosdebye('Fe',300.0)
    x.set_dyninfo_vdosdebye('Al',200.0)
    print( x.create_ncmat() )
    #x.load('dcutoff=0.03').info.dump()
    #x.load('dcutoff=0.03').info.dump()
    x.load('dcutoff=0.2')
    #x.load('dcutoff=0.2').info.dump()
    #mat = x.load()
    #mat.info.dump()
    #unit test: ase bulk('Cu', cubic=True)


    x = NCHelper_SimpleCrystalMaterial()
    x.set_cellsg_cubic( 4.0, spacegroup=225 )
    x.set_atompos( [ ('Al',   0, 0.5, 0.5),
                     ('Al',   0,   0,   0),
                     ('Al', 0.5, 0.5,   0),
                     ('Al', 0.5,   0, 0.5) ] )
    #x.set_dyninfo_freegas('Al')
#    x.set_dyninfo_vdosdebye('Fe',300.0)
    x.set_mean_squared_displacement('Al',0.04)

    x.info.dump()
    #x.xsect(0.0)
    #raise SystemExit
    import PyAna#proper FPE disable
    import numpy as np
    import matplotlib.pyplot as plt
    wls = np.linspace(0.0,10.0,500)
    plt.plot( wls, x.xsect(wl=wls) )
    plt.plot( wls, x.xsect(wl=wls,absn=False) )
    plt.show()

    for dcutoff in [0.1,0.8,2.0]:
        x.set_dcutoff(dcutoff)
        plt.plot( wls, x.xsect(wl=wls) )
    plt.show()

    x.set_dcutoff( None )
    x.set_mean_squared_displacement('Al',0.2)
    y = x.clone()
    y.set_mean_squared_displacement('Al',0.01)
    plt.plot( wls, x.xsect(wl=wls) )
    plt.plot( wls, y.xsect(wl=wls), alpha=0.5 )
    plt.show()

if __name__=='__main__':
    _main()
