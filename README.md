## cmpx_spfunc
Library for computing various special functions (such as Gamma, Hypergeometric, Bessel, 
incomplete Beta, etc.) with complex arguments. While other libraries only accept a few 
arguments to be complex (or none), `cmpx_spfunc` provides functions whose arguments are
all defined in terms of complex values (this sometimes requires analytic continuation).

#### Install dependencies

For the compilation, the code relies on [meson](https://mesonbuild.com/)
and [ninja](https://ninja-build.org/). To install meson and ninja, just run: 
```bash
pip install -r requirements.txt
```

In addition, `HpT-MON` relies on a third party C++ library
[complex_bessel](https://blog.joey-dumont.ca/complex_bessel/).

#### Compile & run the code

Thanks to meson, compiling the code is straightforward:
```bash
meson setup builddir
cd builddir
meson compile
```

This will generate an executable called `test_spfunc`. Then, to run the code,
just type the followinc command:
```bash
./test_spfunc
```

Every time changes are made, the code can be re-compiled by just running `meson compile`
inside the `builddir` directory.


#### System-wide installation

Finally, in case one wants to install the header files and library system-wide, this
can be done by running the following:
```bash
meson install
```
This, by default, will install the header files in `/{prefix}/include/higgs-fo` and
add `cmp_spfunc.pc` to the PKG path.
