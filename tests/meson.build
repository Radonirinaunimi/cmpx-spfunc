project(
  'testcmpfunc', 'cpp',
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

include_directories('../cmpfunc_cpp/include')
inc = include_directories(['../cmpfunc_cpp/include'])

install_headers('catch.hpp')

sources = [
  '../cmpfunc_cpp/src/complex_operators.cpp',
  '../cmpfunc_cpp/src/complex_functions.cpp',
]

testlib = static_library(
  'testfunc',
  sources,
  include_directories : inc,
  dependencies : [glib_dep, cmpx_bes],
  install: false
)

dyexe = executable(
  'unittest', './unittests.cpp',
  include_directories : inc,
  dependencies : [glib_dep],
  link_with : testlib,
  c_args : extra_args
)
