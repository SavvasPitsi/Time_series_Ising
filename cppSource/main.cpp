#include "matplotlibcpp.h"
#include <iostream>
#include <thread>

namespace plt = matplotlibcpp;

void plot1(std::string s)
{   
    // plt::backend("TkAgg");

    //plt::detail::_interpreter::kill();
    
    
}

void plot_fig(std::vector<double>& v,int figNum = 1, std::string s="Title")
{
    plt::figure(figNum);
    plt::clf();
    plt::plot(v);
    plt::title(s);
    plt::show(false);
    plt::pause(0.1);
}

int main() {
    // plt::backend("TkAgg"); // not necessary after PyQt6 install and backend fix
  
    
    double A[] = {1, 3, 2, 4, -5};
    double B[] = {-1, 3, -2, 6, 5};
    
    int n1 = sizeof(A) / sizeof(A[0]);
    int n2 = sizeof(B) / sizeof(B[0]);

    std::vector<double> V1(A,A+n1);
    std::vector<double> V2(B,B+n2);

    plot_fig(V1,2, "Plot #1");
    plot_fig(V2);
    
    // plt::plot();
    plt::show();
    // int x;
    // std::cin >> x;

    // std::cout << "\nThis is figure " << plt::figure() << "\n" << std::endl;
    // std::thread t2(plot1, "This is NOT it");
    // t2.join();
    // plt::show();
    // plt::detail::_interpreter::kill();
    // Py_Finalize();
}
