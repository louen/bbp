#include "CoreCpp/CoreMacros.hpp"
#include "CoreCpp/CoreStrings.hpp"

// Naive implementation of the BBP pi estimation algorithm
double bbp_naive(uint N){
    double result = 0.0;
    double div = 1.0;
    for (uint k = 0; k <= N; ++k) {
        const double term = 
             ((4.0 / (k * 8 + 1)) 
            - (2.0 / (k * 8 + 4))
            - (1.0 / (k * 8 + 5))
            - (1.0 / (k * 8 + 6)));
        
        result += div * term;
        //printf("%d \t %.20e \t %.20e \t %.20e \t %.20e \t %.20e \n", k, result, div, term, div*term, result - M_PI);
        printf("%d \t %.20e \t %.20e \n", k, result, result - M_PI);
        div *= 1.0/16.0;
    }
    return result;
}



int main() {
    bbp_naive(10);
}
