from setuptools import setup
import os

__version__ = "0.0.1"

if not os.path.isfile('cpp/build/parinterp.cpython-310-x86_64-linux-gnu.so'):
    raise Exception("Failed to build a C++ shared object")

setup(
    name="parinterp",
    version=__version__,
    author="Alexey Belkov",
    url="https://github.com/alexeybelkov/parinterp/tree/rework",
    description="Parallel 2D linear interpolation library",
    long_description="",
    py_modules=['cpp/build/parinterp'],
    zip_safe=False,
    python_requires=">=3.10",
)
