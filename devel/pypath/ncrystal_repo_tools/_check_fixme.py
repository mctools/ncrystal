
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
    search_str = 'f' + 'i' + 'x' + 'm' + 'e'
    from .srciter import all_files_iter
    from .dirs import reporoot
    hits = []
    for f in all_files_iter():
        n = 0
        try:
            content = f.read_text()
        except UnicodeDecodeError as e:
            #FIXME: check if file is committed or not?
            raise RuntimeError(f'Binary file: {f}') from e
        print(f)
        for line in content.splitlines():
            if search_str in line.lower():
                n+=1
        if n:
            hits.append( ( n, str(f.relative_to(reporoot)) ) )
    if hits:
        ntot = sum( n for n,f in hits )
        hits.append( (ntot, 'TOTAL' ) )
        hits.sort()
        wn = max( len(str(n)) for n,f in hits )
        for n,f in hits:
            print( ' %s %s'%( str(n).rjust(wn), f ) )
        print()
        print( f"ERROR: A total of {ntot} {search_str}'s found!")
        raise SystemExit(1)


if __name__=='__main__':
    main()
