import os
from distutils.core import setup
from distutils.extension import Extension
from Cython import __version__ as cython_version
from Cython.Build import cythonize
from Cython.Distutils import build_ext

bamtools_include = [os.path.abspath('../bamtools/include'),]
bamtools_lib = [os.path.abspath('../bamtools/lib'),]
source_pattern = 'rmatspipeline/%s.pyx'

# Prefer static BamTools if available; fall back to dynamic with rpath
bamtools_libdir = bamtools_lib[0]
bamtools_a = os.path.join(bamtools_libdir, 'libbamtools.a')
use_static_bamtools = os.path.exists(bamtools_a)

# Common kwargs for the extension
ext_kwargs = dict(
    name='rmats.rmatspipeline',
    sources=[source_pattern % "rmatspipeline"],
    include_dirs=bamtools_include,
    library_dirs=bamtools_lib,
    language="c++",
    extra_compile_args=[
        '-O3', '-funroll-loops',
        '-std=c++11', '-fopenmp',
        '-D__STDC_CONSTANT_MACROS',
        '-D__STDC_LIMIT_MACROS', '-w',
    ],
)

# Base libraries (without bamtools)
base_libs = ['m', 'stdc++', 'z']
extra_link_args = ['-fopenmp']

if use_static_bamtools:
    # Link BamTools statically by pulling in the archive directly
    ext_kwargs['extra_objects'] = [bamtools_a]
    ext_kwargs['libraries'] = base_libs
else:
    # Fall back to dynamic bamtools and embed its location for runtime
    ext_kwargs['libraries'] = base_libs + ['bamtools']
    extra_link_args += ['-Wl,-rpath,' + bamtools_libdir]

ext_kwargs['extra_link_args'] = extra_link_args

asevent_ext = [
    Extension(**ext_kwargs)
]

# https://cython.readthedocs.io/en/latest/src/userguide/migrating_to_cy30.html#exception-values-and-noexcept
compiler_directives = dict()
if cython_version.startswith('3'):
    compiler_directives['legacy_implicit_noexcept'] = True

setup(
    name = 'rmats.rmatspipeline',
    ext_modules = cythonize(asevent_ext, compiler_directives=compiler_directives),
    cmdclass = {'build_ext': build_ext},
)
