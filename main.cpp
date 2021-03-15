#include "CoreCpp/CoreMacros.hpp"
#include "CoreCpp/CoreStrings.hpp"

// Naive implementation of the BBP pi estimator up to N terms
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
        div *= 1.0/16.0;
    }
    return result;
}

// Modular exponentiation : computes base ^ exp % mod 
// Fast computation by square-and-multiply method
//https://en.wikipedia.org/wiki/Modular_exponentiation 
uint modular_exp(uint base, uint exp, const uint mod) {
    CORE_ASSERT(mod > 0, "Modulus must be positive");
    if (mod == 1) {
        return 0;
    }
    uint result = 1;
    base = base % mod;
    for (;exp > 0; exp = exp >> 1) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
    }
    return result;
}

template <uint C>
float bbp_sum(uint n) {
    // Warning : our n in code is actually n-1 in the formulae
    // (so that the 0th digit is 3 etc.)

    // finite sum [0.. n]
    float left = 0.f;
    for (uint k = 0; k < n; ++k) {
        const uint den = 8*k + C;
        const uint num = modular_exp(16, n-1-k, den);
        left += (float(num) / float(den)); 
    }
    
    // partial infinite sum [ n+1 .. ]
    const uint extra_iters = 5;
    float div = 16.0;
    float right = 0.f;
        for (uint k = n ; k < n + extra_iters; ++k) {
        const uint den = 8*k + C;
        right += 1.f / (div * den);
        div *= 16.0;
    }
    return left+right;
}

uint bbp(uint n) {
    // Compute bbp
    const float bbp_result = 4.f * bbp_sum<1>(n) - 2.f * bbp_sum<4>(n) - bbp_sum<5>(n) - bbp_sum<6>(n);
    
    // extract nth digit
    const float scaled_result = bbp_result - (int)(bbp_result) + 1.f;
    return (uint)( 16.f * (scaled_result - (int)(scaled_result)));
}

int main() {
    //bbp_naive(10);

    // Hexadecimal pi for ref
    printf("3.243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89452821e638d01377be5466cf34e90c6cc0ac29b7c97c50dd3f84d5b5b54709179216d5d98979fb1bd1310ba698dfb5ac2ffd72dbd01adfb7b8e1afed6a267e96ba7c9045f12c7f9924a19947b3916cf70801f2e2858efc16636920d871574e69a458fea3f4933d7e0d95748f728eb658718bcd5882154aee7b54a41dc25a59b59c30d5392af26013c5d1b023286085f0ca417918b8db38ef8e79dcb0603a180e6c9e0e8bb01e8a3ed71577c1bd314b2778af2fda55605c60e65525f3aa55ab945748986263e8144055ca396a2aab10b6b4cc5c341141e8cea15486af7c72e993b3ee1411636fbc2a2ba9c55d741831f6ce5c3e169b87931eafd6ba336c24cf5c7a325381289586773b8f48986b4bb9afc4bfe81b6628219361d809ccfb21a991487cac605dec8032ef845d5de98575b1dc262302eb651b8823893e81d396acc50f6d6ff383f442392e0b4482a484200469c8f04a9e1f9b5e21c66842f6e96c9a670c9c61abd388f06a51a0d2d8542f68960fa728ab5133a36eef0b6c137a3be4ba3bf0507efb2a98a1f1651d39af017666ca593e82430e888cee8619456f9fb47d84a5c33b8b5ebee06f75d885c12073401a449f56c16aa64ed3aa62363f77061bfedf72429b023d37d0d724d00a1248db0fead3\n");
    for (uint i = 0; i <= 1000; ++i) {
        uint v = bbp(i);
        printf("%x",v);
        if (i ==0) { printf(".");}
    }
    printf("\n");
}
