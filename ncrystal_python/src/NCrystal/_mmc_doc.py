
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

class DocHelper:
    def __init__(self, toptitle, *, is_wiki ):
        self.is_wiki = is_wiki

        self.o = []
        self.add_title(toptitle,lvl=1)

    def add_title(self,t,*,lvl=2):
        self.add_empty()
        if self.is_wiki:
            assert lvl in (1,2,3)
            self.o.append('%s %s'%('#'*lvl,t))
        else:
            self.o.append(t)
            self.o.append(('=' if lvl==1 else '-')*len(t))

        if self.is_wiki and lvl==1:
            self.o.append('<!-- WARNING: PAGE AUTOMATICALLY GENERATED'
                          ' DO NOT EDIT DIRECTLY -->')
        self.add_empty()

    def add_wrap(self,txt,txt2=None):
        si = '' if txt2 is None else ' '*len(txt)
        txt += txt2 or ''
        if self.is_wiki:
            self.add_line(txt)
        else:
            from textwrap import fill as twf
            self.o += twf( txt, width=80, subsequent_indent=si ).splitlines()

    def add_line(self,line):
        self.o.append(line)

    def add_empty( self ):
        if self.o and self.o[-1]:
            self.o.append('')

    def __iadd__(self,lines):
        self.o += lines

    def lines( self ):
        return self.o

    def __str__(self):
        s = '\n'.join(self.o)
        if not s.endswith('\n'):
            s += '\n'
        return s

def check_keys(d,*keys):
    assert set(d.keys())==set(keys)

def gendoc_geom( **kwargs ):
    """Generate documentation page for MiniMC geometries."""
    from NCrystalDev._common import json_query_cpplayer
    data = json_query_cpplayer(['mmc','cfgdoc','geom'])
    doc = DocHelper('MiniMC geometry cfg-strings', **kwargs)
    check_keys( data,'intro_text','geom_list' ) #fixme: 'examples'!!!
    doc.add_wrap(data['intro_text'])
    doc.add_title('Specific geometries available')
    doc.add_wrap('The available geometries are:')
    nmax = max(len(geom['name']) for geom in data['geom_list'])
    for geom in sorted( data['geom_list'],
                        key = lambda e : e['sort_key'] ):
        check_keys(geom,'name','descr','params','sort_key')
        doc.add_empty()
        if doc.is_wiki:
            title=('* `"%s"`:'%geom['name'])
            titlesp = '    '
            doc.add_line(title)
            doc.add_empty()
            doc.add_wrap('    ' if doc.is_wiki else '  ',geom['descr'] )
        else:
            title=('  "%s": '%geom['name']).rjust(nmax+6)
            titlesp = ' '*len(title)
            doc.add_wrap( title,geom['descr'] )
        doc.add_empty()
        doc.add_line(f'{titlesp}Parameters:')
        if doc.is_wiki:
            titlesp = '    *'
        else:
            titlesp += ' '
        for pname, pdefval, pdescr in geom['params']:
            assert pdefval is not None
            a,b = f'{titlesp} {pname} (default {pdefval}) : ', pdescr
            if doc.is_wiki:
                doc.add_line(a+b)
            else:
                doc.add_wrap(a,b)
    return doc


def gendoc_src( **kwargs ):
    """Generate documentation page for MiniMC sources."""
    from NCrystalDev._common import json_query_cpplayer
    data = json_query_cpplayer(['mmc','cfgdoc','src'])
    doc = DocHelper('MiniMC source cfg-strings', **kwargs)
    check_keys( data,
                'intro_text',
                'src_list',
                'examples',
                'commonpars_energy_descr',
                'commonpars_energy_examples',
                'commonpars_descr',
                'commonpars_list' )
    doc.add_wrap(data['intro_text'])

    def print_example_list( list_ex_descr ):
        for ex,descr in list_ex_descr:
            if doc.is_wiki:
                doc.add_line( '* `"%s"`'%ex )
                doc.add_empty()
                doc.add_line(f'    {descr}')
            else:
                doc.add_line('  "%s":'%ex)
                doc.add_wrap('      ',descr)
            doc.add_empty()

    doc.add_title('Quick usage examples')
    print_example_list(data['examples'])

    doc.add_title('Parameters common to all sources')
    doc.add_wrap(data['commonpars_descr'])
    doc.add_empty()
    assert data['commonpars_energy_descr'].endswith('these examples:')
    doc.add_wrap(data['commonpars_energy_descr'])

    doc.add_empty()
    print_example_list(data['commonpars_energy_examples'])

    doc.add_title('Specific sources available')
    doc.add_wrap('The available sources are:')

    nmax = max(len(src['name']) for src in data['src_list'])
    lcommon = data['commonpars_list']
    lcommon = ', '.join(lcommon[:-1])+f', and {lcommon[-1]}'

    for src in sorted( data['src_list'],
                       key = lambda e : e['sort_key'] ):
        check_keys(src,'name','descr','specific_params','sort_key')
        doc.add_empty()
        if doc.is_wiki:
            title=('* `"%s"`:'%src['name'])
            titlesp = '    '
            doc.add_line(title)
            doc.add_empty()
            doc.add_wrap( '    ' if doc.is_wiki else '  ',src['descr'] )
        else:
            title=('  "%s": '%src['name']).rjust(nmax+6)
            titlesp = ' '*len(title)
            doc.add_wrap( title,src['descr'] )
        doc.add_empty()
        doc.add_line(f'{titlesp}Parameters:')
        doc.add_empty()
        if doc.is_wiki:
            titlesp = '    *'
        else:
            titlesp += ' '
        doc.add_line(f'{titlesp} {lcommon} : Explained above.')
        for pname, pdefval, pdescr in src['specific_params']:
            assert pdefval is not None
            a,b = f'{titlesp} {pname} (default {pdefval}) : ', pdescr
            if doc.is_wiki:
                doc.add_line(a+b)
            else:
                doc.add_wrap(a,b )
    return doc

if __name__=='__main__':
    import sys
    if 'geom' in sys.argv[1:]:
        fct = gendoc_geom
    elif 'src' in sys.argv[1:]:
        fct = gendoc_src
    else:
        assert False
    print(str(fct(is_wiki = 'wiki' in sys.argv[1:])),end='')
