#define _USE_MATH_DEFINES

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <vector>
#include <cmath>
#include <complex>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)



std::vector<double> sin_signal(double frequency, int t_start, int t_end, int num_samples)
{
    std::vector<double> return_table;
    double time = t_end - t_start;
    double sample_distance = time / num_samples;
    for (int i = 1; i < num_samples; ++i)
    {
        return_table.push_back(sin(2 * M_PI * frequency * (t_start + (i * sample_distance))));
    }
    return return_table;
}

std::vector<double> cos_signal(double frequency, int t_start, int t_end, int num_samples)
{
    std::vector<double> return_table;
    double time = t_end - t_start;
    double sample_distance = time / num_samples;
    for (int i = 1; i < num_samples; ++i)
    {
        return_table.push_back(cos(2 * M_PI * frequency * (t_start + (i * sample_distance))));
    }
    return return_table;
}

std::vector<double> square_signal(double frequency, int t_start, int t_end, int num_samples)
{
    std::vector<double> return_table;
    double time = t_end - t_start;
    double sample_distance = time / num_samples;
    for (int i = 1; i < num_samples; ++i)
    {
        if (fmod((t_start + i * sample_distance), (1 / frequency)) / (1 / frequency) > (1 / frequency) / 2) {
            return_table.push_back(1);
        }
        else return_table.push_back(0);
    }
    return return_table;
}

std::vector<double> sawtooth_signal(double frequency, int t_start, int t_end, int num_samples)
{
    std::vector<double> return_table;
    double time = t_end - t_start;
    double sample_distance = time / num_samples;
    for (int i = 1; i < num_samples; ++i)
    {

        return_table.push_back( 0.25 * fmod((t_start + i * sample_distance), (1 / frequency)) / (1 / frequency));

    }
    return return_table;
}

std::vector<std::complex<double>> DFT(std::vector<std::complex<double>> x)
{
    std::vector<std::complex<double>> y;
    std::complex<double> i = (0, 1);
    int N = x.size();
    for (int k = 0; k < N; k++)
    {
        std::complex<double> i(0, 1);
        std::complex<double> temp(0, 0);
        for (int n = 0; n < N; n++)
        {
            temp += x[n] * exp(- i * 2.0 * M_PI * (double(k) / N) * double(n));
        }
        y.push_back(temp);
    }
    return y;
}

std::vector<std::complex<double>> IDFT(std::vector<std::complex<double>> x)
{
    std::vector<std::complex<double>> y;
    double N = x.size();
    for (int n = 0; n < N; n++)
    {
        std::complex<double> i(0, 1);
        std::complex<double> temp(0, 0);
        for (int k = 0; k < N; k++)
        {
            temp += x[k] * (double(1) / N) * exp( i * 2.0 * M_PI * (double(k) / N) * double(n));
        }
        y.push_back(temp);
    }
    return y;
}

std::vector<double> apply_filter(const std::vector<double>& signal, const std::vector<double>& kernel) {
    int n = signal.size();
    int k = kernel.size();
    int pad = k / 2;

    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int j = 0; j < k; ++j) {
            int idx = i + j - pad;
            if (idx >= 0 && idx < n)
                acc += signal[idx] * kernel[j];
        }
        result[i] = acc;
    }

    return result;
}

std::vector<std::vector<double>> apply_filter_2D(
    const std::vector<std::vector<double>>& input,
    const std::vector<std::vector<double>>& kernel)
{
    int rows = input.size();
    int cols = input[0].size();
    int krows = kernel.size();
    int kcols = kernel[0].size();

    int k_center_y = krows / 2;
    int k_center_x = kcols / 2;

    std::vector<std::vector<double>> output(rows, std::vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double acc = 0.0;

            for (int m = 0; m < krows; ++m) {
                for (int n = 0; n < kcols; ++n) {
                    int y = i + m - k_center_y;
                    int x = j + n - k_center_x;

                    if (y >= 0 && y < rows && x >= 0 && x < cols) {
                        acc += input[y][x] * kernel[m][n];
                    }
                }
            }

            output[i][j] = acc;
        }
    }

    return output;
}

std::vector<std::complex<double>> remove_low_f(const std::vector<std::complex<double>>& signal,
int num_to_remove)
{
std::vector<std::complex<double>> spectrum = DFT(signal);

int N = spectrum.size();
for (int i = 0; i < num_to_remove && i < N; ++i) {
    spectrum[i] = 0.0;
}
std::vector<std::complex<double>> filtered = IDFT(spectrum);
return filtered;
}

using namespace matplot;

void plot_signal(const std::vector<double>& y) {
    using namespace matplot;

    // Tworzymy wektor x – osie czasu lub indeksów
    std::vector<double> x(y.size());
    std::iota(x.begin(), x.end(), 0);  // x = [0, 1, 2, ..., N-1]

    plot(x, y);
    title("Wykres sygnału");
    xlabel("Próbki");
    ylabel("Amplituda");
    grid(on);
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
    
    m.def("DFT", &DFT, R"pbdoc(
        Transformata dyskretna
        
    )pbdoc");
    
    m.def("IDFT", &IDFT, R"pbdoc(
        odwrootnoœæ transformaty
    
    )pbdoc");
    
    m.def("sin_signal", &sin_signal, R"pbdoc(
            sygna³ sinus
        
    )pbdoc");
    m.def("cos_signal", &cos_signal, R"pbdoc(
            sygna³ cosinus
        
    )pbdoc");

    m.def("square_signal", &square_signal, R"pbdoc(
        sygna³ prostok¹tny
        
    )pbdoc");

    m.def("sawtooth_signal", &sawtooth_signal, R"pbdoc(
        sygna³ pi³ozêbny
        
    )pbdoc");
    
    m.def("plot_signal", &plot_signal,  R"pbdoc(
        Pokazuje wykres
        
    )pbdoc");
    
    m.def("apply_filter", &apply_filter, R"pbdoc(
            nakladanie filtra 1D
        
    )pbdoc");

    m.def("apply_filter_2D", &apply_filter_2D, R"pbdoc(
            nakladanie filtra 2D na macierz 2D
        
    )pbdoc");

    m.def("remove_low_frequencies_using_DFT", &remove_low_frequencies_using_DFT,
"Usuwa niskie częstotliwości  z widma sygnału",
        
pybind11::arg("signal"), pybind11::arg("num_to_remove"));
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(MACRO_STRINGIFY);
#else
    m.attr("__version__") = "dev";
#endif
}
