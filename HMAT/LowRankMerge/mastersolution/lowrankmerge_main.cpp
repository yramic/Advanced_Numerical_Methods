/**
 * @file lowrankmerge_main.cpp
 * @brief NPDE homework XXX MAIN FILE
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "lowrankmerge.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>

int main(int argc, char **argv) {

    auto start = std::chrono::high_resolution_clock::now();
    // Tabulate approximation errors of low_rank_merge()
    std::cout << "Problem 2-6.c" << std::endl;
    for(int p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,double> errs = LowRankMerge::test_low_rank_merge(n);
        std::cout << "n = " << n << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Max = "
                  << std::scientific << std::setprecision(3)
                  << errs.second << std::endl;
    }

    // Tabulate approximation errors and resulting ranks of adap_rank_merge()
    std::cout << "Problem 2-6.e"       << std::endl;
    std::cout << "Fixed rtol = 0.0001" << std::endl;
    double rtol = 0.0001;
    for(int p=3; p<=12; ++p) {
        unsigned n = std::pow(2,p);
        std::pair<double,size_t> errs = LowRankMerge::test_adap_rank_merge(n, rtol);
        std::cout << "n = " << n << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
    std::cout << "Problem 2-6.e" << std::endl;
    std::cout << "Fixed n = 500" << std::endl;
    unsigned n = 500;
    for(int p=1; p<=8; ++p) {
        double rtol = std::pow(10,-p);
        std::pair<double,size_t> errs = LowRankMerge::test_adap_rank_merge(n, rtol);
        std::cout << "rtol = " << rtol << std::setw(15)
                  << "Frobenius = "
                  << std::scientific << std::setprecision(3)
                  << errs.first  << std::setw(10)
                  << "Rank = "
                  << errs.second << std::endl;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout<<"Runtime: "<<duration.count()<<" ms"<<std::endl;
    return 0;
}