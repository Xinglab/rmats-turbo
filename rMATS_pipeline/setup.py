import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

bamtools_include = [os.path.abspath('../bamtools/include'),]
bamtools_lib = [os.path.abspath('../bamtools/lib'),]
source_pattern = 'rmatspipeline/%s.pyx'

os.environ['CC'] = 'g++'
os.environ['CXX'] = 'g++'
# os.environ['CXXFLAGS'] = '-Wl,-rpath,%s' % (bamtools_lib[0])
# os.environ['LD_LIBRARY_PATH'] = '%s' % (bamtools_lib[0])

asevent_ext = [
    Extension('rmats.rmatspipeline', sources=[source_pattern % "rmatspipeline"],
              include_dirs=bamtools_include,
              libraries=['m','stdc++','bamtools','z'],
              library_dirs=bamtools_lib,
              extra_compile_args = ['-O3', '-funroll-loops',
                                    '-std=c++0x', '-fopenmp',
                                    '-D__STDC_CONSTANT_MACROS',
                                    '-D__STDC_LIMIT_MACROS', '-w',
                                    '-Wl,-static',],
              extra_link_args=['-lbamtools', '-lm', '-std=c++0x',
                               '-lz', '-fopenmp',],
              language="c++",
              )
    ]

setup(
    name = 'rmats.rmatspipeline',
    ext_modules = cythonize(asevent_ext),
    cmdclass = {'build_ext': build_ext},
)
