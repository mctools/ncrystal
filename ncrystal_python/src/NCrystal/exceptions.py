
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

"""Module defining all the exception and warning types emitted by NCrystal"""

#NB: We also put all of these in __all__ in the main __init__.py:

__all__ = [ 'NCrystalUserWarning',
            'NCException',
            'NCFileNotFound',
            'NCDataLoadError',
            'NCMissingInfo',
            'NCCalcError',
            'NCLogicError',
            'NCBadInput',
            'nc_assert' ]

class NCrystalUserWarning( UserWarning ):
    """UserWarning's emitted from NCrystal code"""
    def __init__(self,*args,**kwargs):
        super(NCrystalUserWarning, self).__init__(*args,**kwargs)

class NCException(RuntimeError):
    """Base class for all exceptions raised by NCrystal code"""
    pass

class NCFileNotFound(NCException):
    pass

class NCDataLoadError(NCException):
    pass

class NCMissingInfo(NCException):
    pass

class NCCalcError(NCException):
    pass

class NCLogicError(NCException):
    pass

class NCBadInput(NCException):
    pass

def nc_assert(b,msg=""):
    """Assertion which throws NCLogicError on failure"""
    if not bool(b):
        raise NCLogicError(msg if msg else 'assertion failed')
