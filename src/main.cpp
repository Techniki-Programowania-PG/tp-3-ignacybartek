#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

vector <double> DFT(vector<double> x)
{
    vector<double> y;
    N = x.sizeof();
    for (int k = 0; k < N; k++)
    {
        double temp=0;
        for (int n = 0; n < N; n++)
        {
            temp+= x[n] * exp()
        }
    }
}

using namespace matplot;

void pokaz() {
    plot({ 1, 2, 3, 4 });
    show();
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
    m.def("pokaz", &pokaz,  R"pbdoc(
        "Pokazuje wykres"

        
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
