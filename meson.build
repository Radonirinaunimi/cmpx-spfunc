project(
  'cmp_spfunc', 'cpp',
  version : '1.0.0',
  default_options : ['cpp_std=c++14', 'warning_level=1'],
  license : 'MIT'
)

cxx = meson.get_compiler('cpp')

if cxx.get_id() == 'clang'
  extra_args = ['-fclang-flag']
else
  extra_args = []
endif

glib_dep = [dependency('gsl')]
cmpx_bes = cxx.find_library('complex_bessel', required: true)

inc = include_directories(
  ['include']
)

subdir('include')
subdir('src')

# create executable for main
higgsexen = executable(
  'test_spfunc', 'test_spfunc.cpp',
  include_directories : inc,
  dependencies : [glib_dep, cmpx_bes],
  link_with : spfunclib,
  c_args : extra_args
)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries : spfunclib,
  version : '1.0.0',
  name : 'cmp_spfunc',
  description : 'A library for special functions with complex arguments.'
)
