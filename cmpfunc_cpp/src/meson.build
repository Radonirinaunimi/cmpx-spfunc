sources = [
  './complex_operators.cpp',
  './complex_functions.cpp',
]

spfunclib = static_library(
  'cmp_spfunc',
  sources,
  include_directories : inc,
  dependencies : [glib_dep, cmpx_bes],
  install: true
)
