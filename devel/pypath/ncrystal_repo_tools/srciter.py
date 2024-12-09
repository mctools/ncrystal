
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

import pathlib
import fnmatch

def all_files_iter( *patterns, root = None ):
    patternset = PatternSet( *expand_patterns( patterns ))
    gitignore = get_main_gitignore()
    from .dirs import reporoot
    for f in _all_files_iter_impl( root or reporoot,
                                   patternset,
                                   gitignore ):
        yield f

def _all_files_iter_impl( currentdir, patternset, gitignore ):
    for p in currentdir.iterdir():
        if p.is_symlink():
            continue#always just ignore
        if gitignore.is_ignored(p):
            #Ignore gitignored files, and do not dive into git-ignore'd subdirs.
            continue

        if p.is_dir():
            #We do not yield or pattern test dirs, but we do descend to look for
            #files:
            if p.name=='.git':
                continue#always ignore
            for f in _all_files_iter_impl( p, patternset, gitignore ):
                yield f
        else:
            if not patternset.accepts( p ):
                continue
            if p.name=='.gitignore':
                continue#always ignore
            yield p

class PatternSet:

    def __init__( self, *pattern_expressions ):
        self.__neg_patterns = []
        self.__pos_patterns = []
        for pexpr in pattern_expressions:
            p = SinglePattern(pexpr)
            if p.is_negated():
                self.__neg_patterns.append( p )
            else:
                self.__pos_patterns.append( p )

    def accepts( self, path ):
        for p in self.__neg_patterns:
            if p.do_match( path ):
                return False
        if self.__pos_patterns:
            #positive patterns are present, must satisfy at least one:
            for p in self.__pos_patterns:
                if p.do_match( path ):
                    return True
            return False
        else:
            #no positive patterns are present, by default this means accept:
            return True

class SinglePattern:
    def __init__( self, pattern ):
        if pattern.startswith('!'):
            self._is_negated = True
            self.__pattern = pattern[1:]
        else:
            self._is_negated = False
            self.__pattern = pattern
        if pattern.startswith('<fileinit>'):
            self._match_on_file_start = True
            self.__pattern = self.__pattern[10:]+'*'
        else:
            self._match_on_file_start = False

    def __str__(self):
        return ( 'SinglePattern('
                 '%s, isneg=%s, onfilestart=%s)'%( repr(self.__pattern),
                                                   self._is_negated,
                                                   self._match_on_file_start) )

    def is_negated( self):
        return self._is_negated

    def do_match( self, path ):
        if self._match_on_file_start:
            if not path.is_file():
                return False
            with path.open('rb') as fh:
                first = fh.read(len(self.__pattern)+10).decode('ascii',
                                                               'replace')
                return do_match_pattern( first, self.__pattern )
        else:
            return do_match_pattern( path, self.__pattern )

special_patterns_db = {
    'cpp':['*.hh','*.cc','*.hpp','*.cpp','*.icc'],
    'c': ['*.h','*.c'],
    'py': ['*.py','<fileinit>#!/usr/bin/env python'],
    'bash': ['*.sh',
             '<fileinit>#!/usr/bin/env bash',
             '<fileinit>#!/bin/bash' ],
    'txt':['*.log','*.md','*.txt','*.cfg',
           'INSTALL','README','FILES','CHANGELOG','VERSION'],
    'cmake':['*.cmake','*/CMakeLists.txt','*cmake*.txt'],
    'toml':['*.toml'],
    'yml':['*.yml','*.yaml'],
    'ncmat':['*.ncmat'],
}

def expand_patterns( patterns ):
    match_patterns = []
    for p in patterns:
        if p.startswith('!'):
            p_special_key = p[1:]
            p_special_negated = True
        else:
            p_special_key = p
            p_special_negated = False
        special_patterns = special_patterns_db.get(p_special_key)
        if special_patterns:
            if p_special_negated:
                match_patterns += ['!%s'%e for e in special_patterns ]
            else:
                match_patterns += special_patterns
        else:
            #Not a special key, just add the pattern:
            match_patterns.append( p )
    return match_patterns

def do_match_pattern( path, pattern ):
    case_insensitive = False
    if pattern.endswith('::i'):
        case_insensitive = True
        pattern = pattern[:-3]
    if hasattr(path,'__fspath__'):
        path = pathlib.Path(path)
        s = str(path.absolute())
        if path.is_dir():
            s += '/'
    else:
        s = path
    if case_insensitive:
        s = s.lower()
        pattern = pattern.lower()
    return fnmatch.fnmatchcase(s,pattern)

class GitIgnoreFile:

    """A rough implementation of .gitignore logic. Partially implements the
    logic as described on https://git-scm.com/docs/gitignore , but this is by no
    means a complete and proper implementation - just good enough for our
    purposes.
    """

    def __init__( self, gitignore_file ):
        self.__patterns = _load_gitignore_patterns(gitignore_file)

    def is_ignored( self, path ):
        s = str(path.absolute())
        is_dir = path.is_dir()
        is_matched_overall = False
        for is_negated,  dirs_only, pattern in self.__patterns:
            is_matched = False
            if dirs_only and not is_dir:
                is_matched = False
            else:
                is_matched = do_match_pattern( s, pattern )
            if is_matched:
                is_matched_overall = not is_negated
        return is_matched_overall

def _load_gitignore_patterns( gitignore_file = None):
    from . import dirs
    if gitignore_file is None:
        gitignore_file = dirs.reporoot.joinpath('.gitignore')
    gitignore_file = gitignore_file.absolute()
    assert gitignore_file.is_file()
    root = gitignore_file.parent

    patterns = []
    for line in gitignore_file.read_text().splitlines():
        if not line or line[0]=='#':
            continue
        while not line.endswith('\\ ') and line[-1].isspace():
            line = line[:-1]
        if not line:
            continue
        if line.startswith('\\!'):
            line = line[1:]
        is_negated = False
        if line[0]=='!':
            is_negated = True
            line = line[1:]
        if not line:
            continue
        line=line.replace('\\#','#').replace('\\ ',' ')
        dirs_only = False
        if line.endswith('/'):
            line = line[:-1]
            dirs_only = True
        if not line:
            continue
        if line.startswith('/'):
            pattern = str(root) + line
        else:
            pattern = str(root) + f'/*{line}'
        patterns.append( ( is_negated,  dirs_only, pattern ) )
    return patterns

_main_gi = [None]
def get_main_gitignore():
    if _main_gi[0] is None:
        from .dirs import reporoot
        _main_gi[0] = GitIgnoreFile( reporoot.joinpath('.gitignore') )
    return _main_gi[0]

