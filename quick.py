def parse_scenario_string():
    """Parses a MiniMC quick simulation scenario string, according to the
    syntax:

    "ENERGY [pencil] [DIVERGENCE] on [THICKNESS] [sphere|slab]"

    Here ENERGY is the monochromatic beam energy

    After parsing, the scenario string is converted into three separate
    cfg-strings for the MiniMC source, geometry, and simulation engine. These
    are then returned in a dictionary like:

    { "srccfg" : "<the source cfg-string>",
      "geomcfg" : "<the geometry cfg-string>",
      "enginecfg" : "<the engine cfg-string>" }
    """

    
