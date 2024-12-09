def main( parser ):
    parser.init( """Launch all or some of the various source code checks
    (defined in <reporoot>/devel/pypath/ncrystal_repo_tools/_check_*.py). The
    default is to run all of them.""" )

    parser.add_argument(
        '-l','--list', action='store_true',
        help="""List all available checks and exit."""
    )
    parser.add_argument(
        'CHECK', nargs='*',
        help="""Run this check (runs all checks if nothing is explicitly
        requested)."""
    )
    args = parser.parse_args()
    from . import check_runner
    if args.list:
        for t in check_runner.get_available_checks_list():
            print(t)
        return
    if not args.CHECK:
        check_runner.run_all_checks()
        return
    done = set()
    for c in args.CHECK:
        if c not in done:
            check_runner.run_check(c)
            done.add(c)
