
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

from ._cliimpl import ( create_ArgumentParser,
                        cli_entry_point,
                        print )
import pathlib

def parseArgs( progname, args, return_parser=False ):

    descr="""
Load input file (for instance an .ncmat file) with NCrystal and verify that the
atom-positions are compatible with the indicated space group. Input must be in a
format supported by NCrystal (e.g. .ncmat), and must represent a crystalline
material with space group and atom position information. And exit code of 0
indicates no issues found with atom positions, exit code of 1 indicates issues,
and exit code of 99 indicates that the material is not crystalline with space
group and atom positions.

The symmetry checking is internally performed with the help of the ASE[1]
spacegroup module, so it might be necessary to install this module first with:

  python3 -mpip install ase

References:
[1]: Ask Hjorth Larsen et al 2017 J. Phys.: Condens. Matter 29 27300
     https://doi.org/10.1088/1361-648X/aa680e
"""
    parser = create_ArgumentParser(prog = progname,
                                   description=descr)
    parser.add_argument("FILE",help="File to load")
    parser.add_argument("--quiet",default=False,action='store_true',
                        help="Will not produce any output unless errors are found")
    parser.add_argument("--wyckoff",default=False,action='store_true',
                        help='Output just Wyckoff positions in easily parsable format (implies --quiet, except for the Wyckoff positions).')
    def_eps = 1e-6
    parser.add_argument("--epsilon",default=def_eps,type=float,
                        help="Level of precision required")

    if return_parser:
        return parser
    args = parser.parse_args(args)
    if args.wyckoff:
        args.quiet = True
    return args

def create_argparser_for_sphinx( progname ):
    return parseArgs(progname,[],return_parser=True)

@cli_entry_point
def main( progname, args ):
    args = parseArgs( progname, args )
    from ._common import get_ncrystal_print_fct
    wyckoff_print = get_ncrystal_print_fct() if args.wyckoff else None
    if args.quiet:
        from ._common import modify_ncrystal_print_fct_ctxmgr
        with modify_ncrystal_print_fct_ctxmgr('block'):
            _main_impl(args,wyckoff_print=wyckoff_print)
    else:
        _main_impl(args)

def _main_impl( args, wyckoff_print = None ):

    if wyckoff_print is None:
        wyckoff_print = print

    try:
        import ase.spacegroup
    except ImportError as e:
        raise RuntimeError('Could not import ase.spacegroup. Depending on your'
                           ' environment, this can be installed with'
                           ' "python3 -mpip install ase" or '
                           '"conda install -c conda-forge ase".') from e
    #Load data:
    from .core import createInfo
    info = createInfo(args.FILE)
    print(f'Loaded info based on: "{args.FILE}"')
    if info.isMultiPhase():
        raise RuntimeError('Not applicable: Material has multiple phases')
    if not info.atominfos:
        raise RuntimeError('Not applicable: Material does not have atom positions')
    if not info.hasStructureInfo() or not info.getStructureInfo()['spacegroup']:
        raise RuntimeError('Not applicable: Material does not have unit cell structure info')
    sg_no = info.getStructureInfo()['spacegroup']

    #Verify filename is not misleading:
    for fnpart in pathlib.Path(args.FILE).stem.split('_'):
        for pattern in ('sg','spacegroup','sgno','spacegr','spacegrp'):
            if fnpart.lower().startswith('sg'):
                _ = fnpart[2:]
                if _ and _.isdigit() and int(_)!=sg_no:
                    raise RuntimeError(f"Filename indicates a different spacegroup number ({int(_)}) than the one specified in the data ({sg_no}).")
    # Analysis utils:
    from ._numpy import (_ensure_numpy, _np )
    _ensure_numpy()
    np = _np
    _cell_offsets = np.asarray(list((a,b,c) for a in (-1,0,1) for b in (-1,0,1) for c in (-1,0,1)),dtype=float)
    class UnitCellPoint:
        def __init__(self,xyz):
            """Unit cell point, all coordinates will be integrally shifted to have values in [0,1)"""
            self._pt = np.asarray(xyz,dtype=float)
            assert len(self._pt)==3
            for i in range(3):
                while xyz[i]<0.0:
                    xyz+=1.0
                while xyz[i]>=1.0:
                    xyz-=1.0
            self._pt = np.asarray(xyz,dtype=float)
            self._mperiodicpts = -(self._pt + _cell_offsets)

        def distSq(self,other):
            """Calculate distance-squared to another unit cell point, taking into account
            the periodicity of the lattice (i.e. points at (0.001,0.001,0.001) and
            (0.999,0.999,0.999) are actually very close to each other)."""
            return np.sum(( self._mperiodicpts + other._pt )**2,axis=1).min()

        def __str__(self):
            return str(self._pt)

        def fmt(self,ndec = 7):
            return (f'(%.{ndec}g, %.{ndec}g, %.{ndec}g)')%tuple(self._pt)

    def matchUpPoints(pts1,pts2):
        """Match up points in the two lists of UnitCellPoint's, finding those of
        closest distance. Returns matches,leftovers where matches is a list of
        (pt1,pt2,dist), and leftovers is a list of unmatched points in case lists
        are of unequal length. The argument lists will be consumed.
        """
        n1 = len(pts1)
        n2 = len(pts2)
        assert n1 > 0 and n2 > 0
        assert isinstance(pts1[0],UnitCellPoint) and isinstance(pts2[0],UnitCellPoint)
        #eps=0.001

        #Create a matrix of distances between pts in the two lists:
        dists = np.empty( shape=(n1,n2) )
        for i1,p1 in enumerate(pts1):
            for i2,p2 in enumerate(pts2):
                dists[i1,i2] = p1.distSq(p2)

        #Keep "popping" the pt pairs with the smallest distances out of the data
        #structures, while recording them in the matches list:
        matches = []
        for i in range(min(n1,n2)):
            idx = np.argmin(dists)
            i2, i1 = idx%len(pts2), idx//len(pts2)
            thedist = np.sqrt(dists[i1,i2])
            #delete entries for these pts and record in list of matches:
            if dists.shape[0]>0:
                dists = np.delete(dists, (i1), axis=0)
            if dists.shape[1]>0:
                dists = np.delete(dists, (i2), axis=1)
            p1 = pts1.pop(i1)
            p2 = pts2.pop(i2)
            matches.append( ( p1, p2, thedist ) )
        assert not ( pts1 and pts2 )#one or both should be empty
        return matches, pts1 or pts2 or []

    def analyse(info, spacegroup,sgdescr):
        """Loop over atoms and verify/analyse. Returns False in case of issues"""
        sg=spacegroup
        errors = False
        wyckoff=[]

        matchlvl = args.epsilon
        nprintdecs = max(7,1+int(np.ceil(-np.log10(matchlvl))))

        distmax = 0.0
        for ai in info.atominfos:
            lbl=ai.atomData.displayLabel()
            print(f"\n==========================> Investigating element: {lbl}")
            print(f"\nThe {len(ai.positions)} positions are generated from the following symmetry-unique points:")
            unique_sites = list(e for e in sg.unique_sites(ai.positions))
            tags = sg.tag_sites(unique_sites + list(ai.positions))

            for idx,us in zip(tags[0:len(unique_sites)],unique_sites):
                wyckoff.append( ( lbl, tuple(us) ) )
                print(f"\nSymmetry-unique point: {us}")
                #equiv_sites = sg.equivalent_sites([us])[0]
                expected_sites = list(UnitCellPoint(e) for e in sg.equivalent_sites([us])[0])
                actual_sites = list(UnitCellPoint(p) for i,p in enumerate(ai.positions) if tags[len(unique_sites)+i]==idx)
                matches,leftovers = matchUpPoints( expected_sites, actual_sites )
                for p1,p2,dist in matches:
                    distmax = max(dist,distmax)
                    problem= bool(dist>matchlvl)
                    errors = errors or problem
                    expected_str = p1.fmt(nprintdecs)
                    actual_str = p2.fmt(nprintdecs)
                    problem_str = ' <-- PROBLEM!!!' if problem else ''
                    expected_str = '' if expected_str==actual_str else ' '+expected_str
                    print(f"    -> Point {actual_str} deviates {dist:g} from expected position{expected_str}{problem_str}")
                if leftovers:
                    errors = True
                    if len(expected_sites)>len(actual_sites):
                        print("\n    ERROR: The following expected symmetry-calculated points were absent from the actual list of points:")
                    else:
                        print("\n    ERROR: The following points were not found among the symmetry-calculated points:")
                    for p in leftovers:
                        print('       -> '+p.fmt(nprintdecs))

        print("\n==========================> Done.\n")
        if not errors:
            print(f"All OK! Atom positions are from these Wyckoff positions for {sgdescr}:\n")
            if args.wyckoff:
                wyckoff_print(f'# Wyckoff positions for {sgdescr}:')
            def fmt( x ):
                return (f'%.{nprintdecs}g')%x
            for lbl,pos in wyckoff:
                wyckoff_print(f'{"  " if not args.wyckoff else ""}{lbl} {fmt(pos[0])} {fmt(pos[1])} {fmt(pos[2])}')
            print(f'\nWorst discrepancy in distances: {distmax:.3g}')
            return True
        else:
            return False
        return True

    sg_setting1 = ase.spacegroup.Spacegroup(sg_no,setting=1)
    print(f"Input has space group: {sg_no} ({sg_setting1.symbol})")

    #Some space-groups support an alternate setting (different point of origin?):
    try:
        sg_setting2 = ase.spacegroup.Spacegroup(sg_no,setting=2)
    except ase.spacegroup.spacegroup.SpacegroupNotFoundError:
        sg_setting2 = None

    assert sg_setting1 is not sg_setting2

    if analyse(info, sg_setting1, sgdescr=f'SG {sg_no} (setting=1)'):
        return # all ok
    problem_msg = ('Problems detected in list of atom positions! Most likely'
                   ' this is a real problem, but you can also check with'
                   f' (enter spacegroup {sg_no}): '
                   'https://www.cryst.ehu.es/cryst/get_wp.html')
    if sg_setting2 is None:
        raise RuntimeError(problem_msg)
    print('\n\nProblems encountered. Re-trying with alternative'
          ' space-group setting (ase.spacegroup.Spacegroup(setting=2))\n\n')
    if analyse(info, sg_setting2, sgdescr=f'SG {sg_no} (setting=2)'):
        return # all ok
    raise RuntimeError(problem_msg)
