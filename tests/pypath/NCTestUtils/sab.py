
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

def extract_scathelper_uids( proc ):
    #Search process for any SABScatter instances and return the SABScatterHelper
    #UID's found in a list.
    if isinstance( proc, dict ):
        #already a summary dict
        s = proc
    elif hasattr(proc,'getSummary'):
        # a scatter process
        s = proc.getSummary(short=False)
    else:
        # a loaded material
        s = proc.scatter.getSummary(short=False)

    uids = []
    if s['name'] == 'ProcComposition':
        for frac, comp in s['specific']['components']:
            uids += extract_scathelper_uids( comp )
    if s['name'] == 'SABScatter':
        uids.append( s['specific']['sabhelper_uid'] )
    return uids
