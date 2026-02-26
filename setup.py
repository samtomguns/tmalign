from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

ext = Pybind11Extension(
    "_tmalign",
    ["tmalign_core.cpp"],
    extra_compile_args=["-O3", "-ffast-math"],
)

setup(
    name="tmalign",
    ext_modules=[ext],
    cmdclass={"build_ext": build_ext},
    py_modules=[],
)
