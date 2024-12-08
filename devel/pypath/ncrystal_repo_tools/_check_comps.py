
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

from .core_components import load_components

def check_comps():
    name2comp = load_components()#This tests against cyclic deps
    def fmt( n ):
        return '%s*'%n if name2comp[n].is_internal else n
    w = max(len(n) for n in name2comp)
    def sort_order( a ):
        return ( len(a[1].deps), a[0] )
    for n,c in sorted(name2comp.items(),key=sort_order):
        print('%s DEPENDS %s'%(fmt(n).ljust(w+1),
                               ' '.join(fmt(e) for e in sorted(c.depnames))))
        #FIXME: Actually test something in each comp. Perhaps test the include
        #statements?!

def main():
    check_comps()

if __name__=='__main__':
    main()
