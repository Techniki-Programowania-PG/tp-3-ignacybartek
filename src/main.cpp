#include <pybind11/pybind11.h>
#include <extern/matplotplusplus>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

void plot(vector <double> y, vector <double> x)
{
    plot(x,y)
    
}
namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";
    m.def("plot", &plot, R"pbdoc(
        ploting 2 numbers

        
    )pbdoc");

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
