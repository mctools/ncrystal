
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

licenseblurb_data = """
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
""".lstrip()

#Fixme: currently fails for a lot of files!

def licenseblurb_data_with_hash_comments():
    assert licenseblurb_data.startswith('///////')
    assert licenseblurb_data.endswith('///////\n')
    out = []
    for line in licenseblurb_data.splitlines():
        #replace / with # but only in the frame
        if line.count('/')>40:
            line = line.replace('/','#')
        else:
            assert line.startswith('//') and line.endswith('//')
            line = '##' + line[2:-2] + '##'
        out.append(line)
    out = '\n'.join(out)+'\n'
    assert out.startswith('#######')
    assert out.endswith('#######\n')
    return out

def do_check( files, allowed_file_hdr_list ):
    for f in files:
        text = f.read_text()
        if not text.strip():
            continue#allow empty files to not have license blurbs
        if not any(text.startswith(e) for e in allowed_file_hdr_list):
            print('bad license blurb:',f)

def main():
    from .srciter import all_files_iter
    lbc = '\n' + licenseblurb_data
    lbh = '\n' + licenseblurb_data_with_hash_comments()

    do_check(all_files_iter('py'),[lbh,'#!/usr/bin/env python3\n'+lbh])
    do_check(all_files_iter('cpp','c'),[lbc,lbc.lstrip()] )

if __name__=='__main__':
    main()
