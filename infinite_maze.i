%module infinite_maze

%{
#include "infinite_maze.h"
%}

%include <stdint.i>
%include <std_pair.i>
%include <std_vector.i>

%template(Point) std::pair<int,int>;
%template(Path) std::vector<std::pair<int,int>>;

%include "infinite_maze.h"