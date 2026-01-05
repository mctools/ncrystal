
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

def load_comps():
    from .core_components import load_components
    name2comp = load_components()#This tests against cyclic deps
    def fmt( n ):
        return '%s*'%n if name2comp[n].is_internal else n
    w = max(len(n) for n in name2comp)
    def sort_order( a ):
        return ( len(a[1].deps), a[0] )
    for n,c in sorted(name2comp.items(),key=sort_order):
        print('%s DEPENDS %s'%(fmt(n).ljust(w+1),
                               ' '.join(fmt(e) for e in sorted(c.depnames))))

    print()
    return name2comp

def check_include_statements( name2comp ):
    from .extract_includes import get_include_staments_from_file

    for c in name2comp.values():
        allowed_h = c.allowed_incs_in_pkgincdir()
        allowed_s = c.allowed_incs_in_pkgsrcdir()
        flist =  [ (f,allowed_h) for f in c.hdrfiles ]
        flist += [ (f,allowed_s) for f in c.srcfiles ]
        flist += [ (f,allowed_s) for f in c.local_hdrs ]
        directly_included_comps = set()
        for f,allowed in flist:
            for i in sorted(get_include_staments_from_file(f)):
                if i not in allowed:
                    p = i.split('/')
                    hint, guess_comp = '', None
                    if ( len(p) in (3,4)
                         and i.startswith('NCrystal/')
                         and p[-2] in name2comp ):
                        guess_comp = p[-2]
                    if guess_comp:
                        hint = ( f'. Perhaps component needs "{guess_comp}"'
                                 f' added in {c.depfile}. You might fix with'
                                 ' the command "ncdevtool fixdeps").' )
                    raise SystemExit(f'ERROR: Forbidden #include "{i}" '
                                     f'in {f}{hint}')
                #include was ok, let us use it for proper dep.txt anaysis:
                p = i.split('/')
                if len(p)==1:
                    #Must be a local header include in srcdir:
                    continue
                if i.startswith('NCrystal/internal/'):
                    assert len(p) == 4
                    directly_included_comps.add( p[2] )
                else:
                    assert i.startswith('NCrystal/')
                    assert len(p) == 3
                    directly_included_comps.add( p[1] )
        current = set(c.direct_depnames)

        #Remove self includes:
        if c.name in directly_included_comps:
            directly_included_comps.remove(c.name)

        if current == directly_included_comps:
            continue

        excessive = current - directly_included_comps
        missing = directly_included_comps - current
        if missing:
            print( f'Missing components in {c.depfile.name}: '
                   + ' '.join(missing) )
        if excessive:
            print( f'Excessive components listed in {c.depfile.name}: '
                   + ' '.join(excessive) )
        if missing or excessive:
            raise SystemExit("ERROR: Based on include statements, the"
                             f" file {c.depfile} has issues reported above.")

def main():
    name2comp = load_comps()
    check_include_statements( name2comp )
