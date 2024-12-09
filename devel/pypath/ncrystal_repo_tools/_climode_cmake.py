
def main( parser ):
    parser.init( 'Launch CMake to build the files in ncrystal_core.' )
    #Fixme many more options here, for now hardcoding below
    parser.add_argument(
        '-b','--build-dir', metavar='DIR',
        help="""CMake build directory to use. By default will work in temporary
        directory which is cleaned up afterwards."""
    )
    parser.add_argument(
        '-i','--install-dir', metavar='DIR',
        help="""CMake install directory to use. By default will work in temporary
        directory which is cleaned up afterwards."""
    )
    parser.add_argument(
        '-f','--force', action='store_true',
        help="""This flag will force any existing build and install dirs to be
        removed before running, instead of ending in an error."""
    )
    args = parser.parse_args()

    needs_tmpdir = (not args.install_dir) or ( not args.build_dir )
    runner_args = dict( force = args.force,
                        mode = 'ctest',
                        cmake_flags = None,
                        build_types = ['rel'],
                        generator = 'single',
                        nprocs_bld = 12,
                        nprocs_ctest = 2 )

    def do_work( blddir, instdir ):
        from .cmake import CMakeRunner
        c = CMakeRunner( blddir = blddir,
                         instdir = instdir,
                         **runner_args )
        c.do_cfg()
        c.do_build()
        c.do_ctest()

    if needs_tmpdir:
        from .util import work_in_tmpdir
        with work_in_tmpdir():
            do_work( args.build_dir or './build',
                     args.install_dir or '/install' )
    else:
        do_work( args.build_dir, args.install_dir )
