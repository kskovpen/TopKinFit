%module kinfit
%{
#include "../include/kinfit.h"
%}

%include "stl.i"
%template(_string_list) std::vector< std::string >;
%include "std_vector.i"
%include "std_string.i"

%template(FloatVector) std::vector<float>;

%include "../include/kinfit.h"


