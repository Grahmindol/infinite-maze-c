from setuptools import setup, Extension

maze = Extension(
    "_infinite_maze",
    sources=["infinite_maze_wrap.cpp"],
    include_dirs=["."],
    define_macros=[("INFINITE_MAZE_IMPLEMENTATION", None)],
    language="c++",
)

setup(
    name="infinite-maze",
    version="0.1.0",
    description="Infinite procedural maze generator",
    author="Grahmindol",
    license="MIT",
    url="https://github.com/Grahmindol/infinite-maze-c",
    py_modules=["infinite_maze"],
    ext_modules=[maze],
)