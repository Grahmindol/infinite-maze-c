from setuptools import setup, Extension

maze = Extension(
    "_infinite_maze",
    sources=["infinite_maze.i"],
    swig_opts=["-c++"],
    include_dirs=["."],
    define_macros=[("INFINITE_MAZE_IMPLEMENTATION", None)],
    language="c++",
)

setup(
    name="infinite-maze",
    version="0.1.0",
    description="Infinite procedural maze generator",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Grahmindol",
    license="MIT",
    url="https://github.com/Grahmindol/infinite-maze-c",
    ext_modules=[maze],
    data_files=[
        ("include", ["infinite_maze.h"]),
    ],
)