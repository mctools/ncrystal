
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

def short_description():
    return 'Visualise C++ components'

def graph_to_dot( graph, fix_size = False ):
    s = ''
    s += 'digraph GG {\n'
    s += '    node [ fontsize = "12" ];\n'

    lbl2id = {}
    for n in graph.nodes:
        n_id = len(lbl2id)
        lbl = n['name']
        lbl2id[ lbl ] = n_id
        display_lbl = n.get('lbl',lbl)
        shape = n.get('shape','octagon')
        color = n.get('color','black')
        fontcolor = n.get('fontcolor',color)
        shapecolor = n.get('shapecolor',color)
        ns = (f'    "node_{n_id}" [ fontcolor="{fontcolor}" '
              f'color="{shapecolor}" label="{display_lbl}" '
              f'shape={shape} ')
        fix_size = False
        if n.get('width') is not None:
            fix_size = True
            ns += 'width=%s '%n.get('width')
        if n.get('height') is not None:
            fix_size = True
            ns += 'height=%s '%n.get('height')
        if fix_size:
            ns += 'fixedsize=true '
        ns += '];\n'
        s += ns

    for a,b,info in graph.connections:
        id_a, id_b = lbl2id[a], lbl2id[b]
        color = (info or {}).get( 'color', 'black' )
        s += f'    "node_{id_a}" -> "node_{id_b}" [ color="{color}" ]\n'

    s += '}\n'
    return s

def render_dot_file( dot_str, fmt='png' ):

    import shutil
    import subprocess
    import pathlib
    from .util import work_in_tmpdir

    cmd_dot = shutil.which('dot')
    cmd_unflatten = shutil.which('unflatten')
    if not cmd_dot:
        raise SystemExit('command "dot" not found (needs graphviz package')

    if not cmd_unflatten:
        raise SystemExit('command "unflatten" not found (needs graphviz package')

    with work_in_tmpdir():
        f1, f2, f3 = './graph1.dot', './graph2.dot','./result'
        pathlib.Path(f1).write_text(dot_str)
        subprocess.run([cmd_unflatten,'-l3','-c7',f1,'-o',f2],check = True )
        f3 = pathlib.Path(f'./graph.{fmt}')
        subprocess.run([cmd_dot,f2,f'-T{fmt}',f'-o{f3}' ],check=True,
                       capture_output = True )
        return pathlib.Path(f3).read_bytes()

def display_image_data( data, fmt ):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from io import BytesIO
    fh = BytesIO(data)
    img = mpimg.imread(fh,fmt)
    plt.imshow(img)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

class Graph:

    def __init__(self):
        self.__nodes = dict()
        self.__connections = dict()

    def add_node( self, nodename, **options ):
        #assert fmt in ['default','disabled']
        assert nodename not in self.__nodes
        self.__nodes[nodename] = options

    def add_connection( self, nodename_mother, nodename_daughter, **options ):
        c = ( nodename_mother, nodename_daughter )
        assert c not in self.__connections
        self.__connections[ c ] = options

    @property
    def nodes( self ):
        for n,o in sorted( self.__nodes.items() ):
            d = { 'name' : n }
            d.update(o)
            yield d

    @property
    def connections( self ):
        return sorted( (a,b,info)
                       for (a,b),info
                       in sorted(self.__connections.items()) )

    def get_as_dot_content( self ):
        return graph_to_dot( self )

def main( parser ):
    parser.init( 'Visualise src component dependencies.' )
    allowed_stats = ('none','sloc','files','headers','sloc_headers')
    allowed_stats_str = '", "'.join(allowed_stats)
    parser.add_argument(
        '-s','--stat', metavar='STAT',
        default = 'none',
        choices=allowed_stats,
        help=('Show stat-count in the visualisation (options:'
              f' "{allowed_stats_str}").'),
    )
    parser.add_argument(
        '-a','--area', action='store_true',
        help="""Use with --stat to make the area of component graphics scale
        with the chosen statistics count."""
    )

    args = parser.parse_args()

    from .core_components import load_components
    import math
    graph = Graph()
    n2c = load_components()
    if args.stat == 'sloc':
        def countfct( c ):
            return c.sloc_count()
    elif args.stat == 'sloc_headers':
        def countfct( c ):
            return c.sloc_count(headers_only = True)
    elif args.stat == 'files':
        def countfct( c ):
            return sum( 1 for f in c.all_file_iter() )
    elif args.stat == 'headers':
        def countfct( c ):
            return len(c.hdrfiles or [])
    else:
        assert args.stat == 'none'
        def countfct( c ):
            return 1

    counts = {}
    for n,c in n2c.items():
        counts[n] = countfct( c )
    if args.area:
        max_count = max( counts.values() )
        min_count = min( counts.values() )
        if max_count and not min_count:
            min_count = min( c for c in counts.values() if c > 0.0 )
        def node_size_fct( c ):
            if max_count == min_count:
                return 0.5
            v = max(min_count,counts[c.name])
            eps = 0.01
            return math.sqrt(eps + (1.0-eps)*( (v-min_count)/(max_count-min_count) ))

    for n,c in n2c.items():
        count = counts[n]
        if args.area:
            node_size = node_size_fct( c )
        else:
            node_size = None
        graph.add_node( n,
                        lbl = f'{n} {count}' if args.stat!='none' else n,
                        color = 'grey' if c.is_internal else 'red',
                        shape = 'circle' if args.area else 'octagon',
                        width = None if node_size is None else 2.0*node_size,
                        height = node_size,
                       )
        for d in c.calc_minimal_deps():
            graph.add_connection(d.name,n, color = 'grey')

    #print( graph.get_as_dot_content())
    fmt = 'png'
    data = render_dot_file( graph.get_as_dot_content(), fmt )
    display_image_data( data, fmt )
    return
#
#    try:
#        import networkx as nx
#    except ImportError:
#        raise SystemExit('ERROR: Missing dependency "networkx"')
#
#    graph = nx.DiGraph()
#    n2c = load_components()
#    for n,c in n2c.items():
#        graph.add_node(n)
#        for d in c.direct_deps:
#            graph.add_edge(d.name,n)
#
#    nx.write_network_text(graph)
#    #return
#
#
#    import matplotlib.pyplot as plt
#    nx.draw(graph,with_labels=True)
#    plt.show()
#    return
#    #nx.draw_random(graph,with_labels=True)
#    #plt.show()
#            #nx.write_network_text(graph)
#    options = {
#        "font_size": 12,
#        "node_size": 3000,
#        "node_color": ((0,0,1,0.2),),
#        "edgecolors": (0,0,0,0.5),
#        "linewidths": 1,
#        "width": 3,
#    }
#    nx.draw_networkx(graph, **options)
#
#    # Set margins for the axes so that nodes aren't clipped
#    ax = plt.gca()
#    ax.margins(0.20)
#    plt.axis("off")
#    plt.show()
#



#    write_network_text(graph, path=None, with_labels=True, sources=None, max_depth=None, ascii_only=False, end='\n', vertical_chains=False)[source]#
