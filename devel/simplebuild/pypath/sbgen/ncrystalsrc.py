
def load_components():
    import sys
    from . import dirs
    e = str(dirs.reporoot/'tests/standalone/pypath')
    if not sys.path or sys.path[0] != e:
        sys.path.insert(0,e)
    import common.ncrystalsrc as mod
    return mod.load_components()
