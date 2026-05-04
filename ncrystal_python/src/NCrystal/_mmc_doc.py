
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

    def list_bullet( self ):
        return '* ' if self.is_wiki else '  * '

    def add_simple_list(self,*entries):
        self.add_empty()
        lb = self.list_bullet()
        for e in entries:
            self.add_wrap(lb,' '.join(e.split()))

    def add_param_list_entry( self, name, descr, *,
                              descr_short = None,
                              defval = None,
                              exvals = None,
                              line_prefix='' ):
        if descr_short is not None and self.is_wiki:
            self.add_line('<details>')
            self.add_line(f'<summary>{name} : {descr_short}</summary>')
            self.add_empty()#required for markdown to work in html tags
            magiceol = '  '#required!
            defval_shown = ( '"%s"'%defval
                             if defval is not None
                             else '(no default value)' )
            self.add_line(f'> **Default value:** {defval_shown}{magiceol}')
            self.add_line(f'> **Description:** {descr}{magiceol}')
            if exvals is not None:
                exvals_shown = '"%s"'%('", "'.join(exvals))
                self.add_line(f'> **Example values:** {exvals_shown}{magiceol}')

            self.add_empty()#required for markdown to work in html tags
            self.add_line('</details>')
            return

        a = f'{line_prefix}{name}'
        if defval is not None:
            a += f' (default "{defval}")'
        a += ' :'
        if self.is_wiki:
            self.add_line(a+' '+descr)
        else:
            self.add_wrap(a+' ',descr)
        if exvals is not None:
            self.add_wrap( ' '*(len(a)+1)+'Example values: ',
                           '"%s"'%('", "'.join(exvals)))

def check_keys(d,*keys):
    assert set(d.keys())==set(keys)

def gendoc_engine( data, **kwargs ):
    doc = DocHelper('MiniMC engine cfg-strings', **kwargs)
    check_keys( data,'intro_text','cfgparams','tallyinfo')
    doc.add_wrap(data['intro_text'])
    d_params = data['cfgparams']
    pnames = sorted(data['cfgparams'].keys())
    doc.add_title('Parameter reference')
    for pn in pnames:
        pd = d_params[pn]
        assert pd['name']==pn and len(pd)==5
        defval = pd['default_value']# '0'
        exvals = pd['example_values']# ['0', '1']
        descr_s = pd['descr_short']
        assert len(descr_s)<=70
        descr_l = pd['descr_long']
        doc.add_param_list_entry( pn, descr_l,
                                  defval=defval,
                                  descr_short = descr_s,
                                  exvals = exvals,
                                  line_prefix = '  ' )
        doc.add_empty()

    doc.add_title('Tally overview')
    hi = data['tallyinfo']['hists']
    h = set(hi.keys())
    doc.add_wrap(' '.join(f"""Currently {len(h)} different quantities are
    available for tallying in histograms as the neutron leaves the geometry
    during a MiniMC simulation. The "tally" parameter is used to select which
    such histograms to produce. Refer also to the parameters "tallybins",
    "tallyref", and "tallybreakdown" which all affects the details of how
    quantities are defined and histogrammed.""".split()))
    doc.add_empty()
    doc.add_line("The available tally quantities are:")
    ll = []
    maxe = 1 if doc.is_wiki else max(len(e) for e in h)
    tt='`' if doc.is_wiki else ''
    for e in sorted(h):
        line = f'{tt}{e.rjust(maxe)}{tt} : {hi[e]["short_descr"]}'
        unit = hi[e]["unit"]
        if unit:
            line += f' ({unit})'
        ll.append(line)
    if doc.is_wiki:
        doc.add_simple_list(*ll)
    else:
        for e in ll:
            doc.add_line(e)
    return doc

def gendoc_geom( data, **kwargs ):
    doc = DocHelper('MiniMC geometry cfg-strings', **kwargs)
    check_keys( data,'intro_text','geom_list','examples' )
    doc.add_wrap(data['intro_text'])
    doc.add_title('Geometry examples')
    doc.add_wrap('Here follows a few examples of how one can set the geometry using the geomcfg string parameter. Note that geometry length values are always to be specified in meters:')
    doc.add_empty()
    _print_example_list(doc,data['examples'])
    doc.add_title('Specific geometries available')
    doc.add_wrap('The full list of available geometries and their parameters:')
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
            doc.add_param_list_entry( pname, pdescr,
                                      defval=pdefval,
                                      line_prefix = f'{titlesp} ' )
    return doc

def _print_example_list( doc, list_ex_descr ):
    for ex, descr in list_ex_descr:
        if doc.is_wiki:
            doc.add_line( '* `"%s"`'%ex )
            doc.add_empty()
            doc.add_line(f'    {descr}')
        else:
            doc.add_line('  "%s":'%ex)
            doc.add_wrap('      ',descr)
        doc.add_empty()


def gendoc_src( data, **kwargs ):
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


    doc.add_title('Quick usage examples')
    _print_example_list(doc,data['examples'])

    doc.add_title('Parameters common to all sources')
    doc.add_wrap(data['commonpars_descr'])
    doc.add_empty()
    assert data['commonpars_energy_descr'].endswith('these examples:')
    doc.add_wrap(data['commonpars_energy_descr'])

    doc.add_empty()
    _print_example_list(doc,data['commonpars_energy_examples'])

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
            doc.add_param_list_entry( pname, pdescr,
                                      defval=pdefval,
                                      line_prefix = f'{titlesp} ' )
    return doc

doc_subjects = ['geom','src','engine','scenario']
def gen_doc_impl( subject, mode ):
    #Notice: Keep docstring of minimc.gen_doc function synchronised with the
    #        implementation here!!!
    from .misc import evaluate_query

    is_wiki = False
    if isinstance(mode,str) and mode.startswith('wiki::'):
        #hidden option
        is_wiki = True
        mode=mode[6:]

    if subject not in doc_subjects:
        from .exceptions import NCBadInput
        raise NCBadInput(f'Unknown subject "{subject}" (must be one of'
                         ' "geom", "src", or "engine")')
    if mode not in ('print','lines','txt','dict'):
        from .exceptions import NCBadInput
        raise NCBadInput(f'Invalid mode "{mode}" (must be one of'
                         ' "print", "lines", "txt" or "dict")')

    if subject == 'scenario':
        if mode == 'dict':
            #Not really supported:
            return {}
        s = gendoc_scenario( is_wiki = is_wiki )
    else:
        data = evaluate_query(['mmc','cfgdoc',subject])
        if mode == 'dict':
            return data
        fctmap = dict( geom = gendoc_geom,
                       src = gendoc_src,
                       engine = gendoc_engine )
        fct = fctmap.get(subject)
        assert fct is not None
        s = fct( data=data, is_wiki=is_wiki )
    if mode=='print':
        from ._common import print as ncprint
        ncprint( str(s), end='' )
    elif mode=='txt':
        return str(s)
    else:
        assert mode=='lines'
        return s.lines()

#Examples, in form of (description, cfgstr, scenariocfg, key_for_test):
_scenariocfg_examples = [
    (
        """Pencil beam of neutrons impinging centrally on a diameter=1mfp (mean
        free path between scatterings) sphere of aluminium. Beam energy is
        chosen to be hopefully interesting for the material, based on Bragg
        threshold and temperature.""",
        "Al_sg225.ncmat",
        "",
        'empty1'
    ),
    (
        """Pencil beam of 2Aa neutrons impinging centrally on a diameter=2mm
        sphere of 300K aluminium.""",
        "Al_sg225.ncmat;temp=300K",
        "2Aa pencil on 2mm sphere",
        'wlpnclonsph'
    ),
    (
        """A zero-divergence beam of 10meV neutrons uniformly illuminating a
        diameter=2mm sphere of 80K beryllium.""",
        "Be_sg194.ncmat;temp=80K",
        "10meV on 2mm sphere",
        'eonsph'
    ),
    (
        """1.8Aa neutrons impinging at right incidence on an infinite slab of
        thickness 10cm filled with humid air.""",
        "gasmix::air/0.9relhumidity",
        "1.8Aa on 10cm slab",
        'wlonslab'
    ),
    (
        """100000 neutrons at a wavelength which is 99% of the Bragg threshold
        of PG, uniformly illuminating a PG filled sphere whose diameter is 2
        times the mean free path length between scatterings.""",
        "C_sg194_pyrolytic_graphite.ncmat",
        "0.99BT on 2mfp 1e5 times",
        'btonmfp'
    ),
]

def gendoc_scenario( **kwargs ):
    doc = DocHelper('MiniMC scenario cfg-strings', **kwargs)
    def pg(s):
        #add paragraph
        doc.add_empty()
        doc.add_wrap(' '.join(s.split()))
    pg("""Under the hood, a MiniMC simulation always needs configuration
          for four separate components:""")
    doc.add_simple_list(
        """A material, defined by a standard NCrystal material cfg-string.""",
        """A source of neutrons, defined by a source cfg-string.""",
        """A simulation geometry, defined by a geometry cfg-string.""",
        """Configuration of the simulation engine itself, defined by an engine
           cfg-string.""")
    pg("""Often the full flexibility offered by these four cfg-strings is not
          needed.  Therefore, as a convenient alternative, one can alternatively
          set up simpler "simulation scenarios", in which configuration of
          source and geometry is automatically handled. No matter what, a
          material cfg-string is always required of course, and an enginecfg can
          still be provided as desired.""")

    doc.add_title('Examples')
    pg("""The following examples show how the combination of two strings are
    enough to setup a MiniMC simulation: one standard NCrystal cfg-string
    defines the material, and one additional string defining the simulation
    scenario.""")

    word_srccfg='[[srccfg|minimc_src]]' if doc.is_wiki else 'srccfg'
    word_geomcfg='[[geomcfg|minimc_geom]]' if doc.is_wiki else 'geomcfg'
    word_enginecfg= ( '[[enginecfg|minimc_engine]]'
                      if doc.is_wiki else 'enginecfg' )
    pg(f"""In each example is also shown the actual underlying {word_geomcfg}
           and {word_srccfg} which are automatically generated based on the
           scenario. An {word_enginecfg} can also be provided separately as
           desired.""")

    bullet = doc.list_bullet()
    bullet_space = ' '*len(bullet)
    v='`' if doc.is_wiki else ''
    doc.add_empty()
    from .minimc import decode_scenario
    for descr, matcfg, scenariostr, key in _scenariocfg_examples:
        dec = decode_scenario( matcfg, scenariostr )
        assert set(dec.keys())==set(['geomcfg','srccfg'])
        gc,sc = dec['geomcfg'],dec['srccfg']
        doc.add_line(bullet      +f'Input:  material={v}"{matcfg}"{v}')
        doc.add_line(bullet_space+f'        scenario={v}"{scenariostr}"{v}')
        doc.add_wrap(bullet_space+'What:   ',' '.join(descr.split()))
        doc.add_line(bullet_space+f'Output: geomcfg={v}"{gc}"{v}')
        doc.add_line(bullet_space+f'        srccfg={v}"{sc}"{v}')

    doc.add_title('Full syntax')
    pg("""The decodeScenario function parses a MiniMC quick simulation scenario
          string, according to the syntax:""")
    syntax = 'ENERGY [pencil] [on [THICKNESS] [sphere|slab]] [COUNT times]'
    doc.add_empty()
    if doc.is_wiki:
        doc.add_line('```')
        doc.add_line(syntax)
        doc.add_line('```')
    else:
        doc.add_line('  '+syntax)

    pg("""The meaning of the various parameters (in uppercase above) and
          keywords (in lowercase above) is explained in the following.""")
    pg("""If given, COUNT is the initial number of neutrons to input into the
          simulation. Otherwise, the default value is 1e6 for isotropic
          materials, and 1e5 for anisotropic materials.""")

    doc.add_simple_list(
        """ENERGY is the monochromatic beam energy like "1.8Aa", "25meV" or
           "0.1eV". Special units "BT" means the bragg threshold of the material
           (or 4.0Aa in case material does not have one), and "kT" means a
           kinetic energy equal to Boltzmann's constant times the material
           temperature. Using either special unit will cause the result to be
           rounded to 6 significant digits.""",
        """"pencil" is an optional keyword related to the beam profile (see
            below).""",
        """THICKNESS is the material thickness like "1mm", "2m", "0.4cm", or
           "2.5mfp". The unit "mfp" corresponds to the mean-free-path length for
           a neutron scattering interaction in the material (rounded to 6
           significant digits). Default THICKNESS is "1mfp".""",
        """The keywords "sphere" or "slab" can be used to select the sample
           geometry (default is a sphere).""",
    )
    pg("""For a spherical geometry, the beam profile will by default be taken to
          be a beam with a uniform circular profile, of the same radius as the
          sphere. However, if the keyword "pencil" is provided, a pencil beam
          hitting the sphere centrally is used instead.""")
    pg("""As a special case, an empty scenario string is interpreted in the same
          way as a scenario string with contents "0.8BT" if the material has a
          Bragg threshold, otherwise it will be "1kT" if it has a
          temperature. If it neither has a Bragg threshold or a temperature, a
          value of "1.8Aa" is used as the ultimate fallback.""")
    pg("""For flexibility and usage from the cmdline, colons (:) and underscores
          (_) can be used as whitespace. Additionally, all repeated whitespace
          (tabs, newlines, etc.) is converted into a single space before
          parsing, and trailing or leading whitespace is trimmed away.""")
    return doc

if __name__=='__main__':
    import sys
    assert len(sys.argv) in (2,3)
    mode = sys.argv[2] if len(sys.argv)==3 else 'print'
    r=gen_doc_impl( subject=sys.argv[1],
                    mode = mode )
    if mode == 'dict':
        from ._common import ncpprint
        ncpprint(r,do_sort=False)
