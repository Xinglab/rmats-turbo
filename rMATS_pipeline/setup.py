import os
from distutils.core import setup
from distutils.extension import Extension
from Cython import __version__ as cython_version
from Cython.Build import cythonize
from Cython.Distutils import build_ext

bamtools_include = [os.path.abspath('../bamtools/include'),]
bamtools_lib = [os.path.abspath('../bamtools/lib'),]
source_pattern = 'rmatspipeline/%s.pyx'

asevent_ext = [
    Extension('rmats.rmatspipeline', sources=[source_pattern % "rmatspipeline"],
              include_dirs=bamtools_include,
              libraries=['m','stdc++','bamtools','z'],
              library_dirs=bamtools_lib,
              extra_compile_args = ['-O3', '-funroll-loops',
                                    '-std=c++11', '-fopenmp',
                                    '-D__STDC_CONSTANT_MACROS',
                                    '-D__STDC_LIMIT_MACROS', '-w',
                                    '-Wl,-static',],
              extra_link_args=['-lbamtools', '-lm', '-std=c++11',
                               '-lz', '-fopenmp',],
              language="c++",
              )
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
