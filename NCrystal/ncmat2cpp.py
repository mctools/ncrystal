
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

def ncmat2cpp( *input_files_or_text_data,
               quiet = False,
               cppfunctionname='registerData',
               compact=False,
               width=140,
               validate=False,
               extra_includes=None,
               outfile = None,
               regfctname=None ):
    """Function which can be used to embed the content of .ncmat files (or actually
    any ASCII/UTF8 excoded text files) directly into a C++ library. It does so
    by reading the files and creating C++ code which keeps the contents of the
    files in static strings, and registers those strings with the NCrystal C++
    API, using the original filename as key. Naturally, the resulting C++ code
    should be stored in a file, and that file must be compiled along with the
    rest of the users C++ code, and the enclosing function must be invoked.

    If outfile is not None, the contents will be stored in that file. In any
    case, the C++ code will be returned as a string.

    if regfctname is None, the C++ function used for registering the in-memory
    file data will be assumed to be
    "NCrystal::registerInMemoryStaticFileData(const std::string&,const char*)".

    FIXME EXPLAIN PARAMETERS (see _cli_ncmat2cpp argparse for descriptions). Be
    sure to point out that despite the name of the module (ncmat2cpp), the
    functionality can actually be used for any UTF-8 encoded text file.

    """
    #NOTE: The above doc-string should be kept in sync with argparse help text
    #in _cli_ncmat2cpp.py.

    #NOTE: ncrystal_ncmat2cpp is a special command-line script, since the
    #NCrystal CMake code needs to be able to invoke _cli_ncmat2cpp.py directly
    #as a standalone script, in a mode where no other of the NCrystal python
    #modules are imported. Therefore, all of the actual code implementing the
    #ncmat2cpp functionality needs to reside in _cli_ncmat2cpp.py, rather than
    #here in ncmat2cpp.py. This is the opposite of how most other command-line
    #scripts are implemented, with the bulk of the implementation residing in
    #the Python API, and the command-line script being a wrapper around it.

    if not input_files_or_text_data:
        from .exceptions import NCBadInput
        raise NCBadInput('No files or text data provided')

    from ._ncmat2cpp_impl import files2cppcode
    return files2cppcode( infiles = input_files_or_text_data,
                          quiet = quiet,
                          outfile = outfile,
                          cppfunctionname = cppfunctionname,
                          compact = compact,
                          width = width,
                          validate = validate,
                          extra_includes = extra_includes,
                          regfctname = regfctname )
