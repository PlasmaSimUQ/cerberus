-- See https://luacheck.readthedocs.io/en/stable/warnings.html
-- 111: Setting undefined global variable (most of our config is globals)
-- 131: unused implicitly defined global variable (some globals are read by the c++ code and may not be used in the lua inputs)

ignore = { "111", "131" }

std = {
  globals = {
    "ref_length",
    "ref_density",
    "ref_mass",
    "ref_temp",
    "lightspeed",
    "beta",
    "skin_depth",
    "Larmor",
    "Debye",
    "verbosity",
    "linear_solver_verbosity",
    "cfl",
    "refine_cutcells",
    "time_integration_scheme",
    "force_dt",
    "refine_boxes",
    "tile_size",
    "plot",
    "states",
    "actions",
    "embedded_boundaries",
    "spairs",
    "len",
    "get_sorted_keys",
    "expand",
    "exists",
    "is_function",
    "is_table",
    "get_ipairs",
    "get_args",
    "get_num_args",
    "preprocess",
    "math",
  },
}

allow_defined_top = true
