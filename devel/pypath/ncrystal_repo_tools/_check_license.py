
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

def do_check( files, allowed_file_hdr_list, *,
              allow_include_guards_if_hdr = False ):

    def ignore_path( frel ):
        #Symlinked NCCFileUtils.* C files are really C++ in their default
        #location, so just skip these files (otherwise we get spurious failures
        #if the symlinks are broken for some reason):
        return frel in [ 'tests/src/app_selfpathc/NCCFileUtils.h',
                         'tests/src/app_selfpathc/NCCFileUtils.c' ]

    def text_ok( text ):
        return any(text.startswith(e) for e in allowed_file_hdr_list)

    def ok_with_incguard( f, text ):
        suffix = f.suffix
        if suffix.endswith('.in') and '.' in f.stem:
            suffix = '.' + f.stem.split('.')[-1]
        if not allow_include_guards_if_hdr or f.suffix not in ('.h','.hh'):
            return False
        lines = text.splitlines()
        if ( len(lines) > 10
             and lines[0].startswith('#ifndef ')
             and lines[1].startswith('#define ') ):
            return text_ok( '\n'.join(lines[2:]) )
        return False

    all_ok = True
    for f in files:
        text = f.read_text()
        if not text.strip():
            continue#allow empty files to not have license blurbs
        if not text_ok(text) and not ok_with_incguard(f,text):
            from .dirs import reporoot
            frel = str(f.relative_to(reporoot)).replace('\\','/')
            if ignore_path( frel ):
                print(f'  Ignoring known false positive: {frel}')
                continue
            print('bad license blurb:',f)
            all_ok = False
    return all_ok

def main():
    from .srciter import all_files_iter
    from .license import ( licenseblurb_data_with_doubledash_comments,
                           licenseblurb_data_with_hash_comments,
                           licenseblurb_data_with_ansic_comments )
    lbcpp = '\n' + licenseblurb_data_with_doubledash_comments()
    lbh = '\n' + licenseblurb_data_with_hash_comments()
    lbansic = '\n' + licenseblurb_data_with_ansic_comments()

    c1 = do_check( all_files_iter('py'),
                   [lbh,'#!/usr/bin/env python3\n'+lbh])
    c2 = do_check( all_files_iter('cpp'),
                   [lbcpp,lbcpp.lstrip()],
                   allow_include_guards_if_hdr = True )
    c3 = do_check( all_files_iter('c'),
                   [lbansic,lbansic.lstrip()],
                   allow_include_guards_if_hdr = True )
    if not (c1 and c2 and c3):
        raise SystemExit('ERROR: license blurb issues detected.')
    print("All ok!")

if __name__=='__main__':
    main()
