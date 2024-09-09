from setuptools import setup
from setuptools import Extension

setup(name='tem_cpp_submodule',
      author='Shahar Sandhaus',
      version='0.0.1',
      packages=["tem_cpp_submodule"],
      ext_modules=[
          Extension('tem_cpp_submodule_native',
                   sources=[
                       'tem_cpp_submodule_native/wrapper.pyx', 
                       'tem_cpp_submodule_native/algs.cpp', 
                       'tem_cpp_submodule_native/seperated_points.cpp'
                   ],
                   language='c++')
        ],
      zip_safe=False)
