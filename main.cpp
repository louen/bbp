#include "CoreCpp/CoreMacros.hpp"
#include "CoreCpp/CoreStrings.hpp"

#include <chrono>
// Naive implementation of the BBP pi estimator up to N terms
double bbp_naive(uint N)
{
    double result = 0.0;
    double div = 1.0;
    for (uint k = 0; k <= N; ++k)
    {
        const double term =
            ((4.0 / (k * 8 + 1)) - (2.0 / (k * 8 + 4)) - (1.0 / (k * 8 + 5)) - (1.0 / (k * 8 + 6)));

        result += div * term;
        div *= 1.0 / 16.0;
    }
    return result;
}

// Modular exponentiation : computes base ^ exp % mod
// Fast computation by square-and-multiply method
//https://en.wikipedia.org/wiki/Modular_exponentiation
uint modular_exp(uint base, uint exp, const uint mod)
{
    CORE_ASSERT(mod > 0, "Modulus must be positive");
    if (mod == 1)
    {
        return 0;
    }
    uint result = 1;
    base = base % mod;
    for (; exp > 0; exp = exp >> 1)
    {
        if (exp % 2 == 1)
        {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
    }
    return result;
}

template <uint C>
float bbp_sum(uint n)
{
    // Warning : our n in code is actually n-1 in the formulae
    // (so that the 0th digit is 3 etc.)

    // finite sum [0.. n]
    float left = 0.f;
    for (uint k = 0; k < n; ++k)
    {
        const uint den = 8 * k + C;
        const uint num = modular_exp(16, n - 1 - k, den);
        left += (float(num) / float(den));
    }

    // partial infinite sum [ n+1 .. ]
    const uint extra_iters = 5;
    float div = 16.0;
    float right = 0.f;
    for (uint k = n; k < n + extra_iters; ++k)
    {
        const uint den = 8 * k + C;
        right += 1.f / (div * den);
        div *= 16.0;
    }
    return left + right;
}

uint bbp(uint n)
{
    // Compute bbp
    const float bbp_result = 4.f * bbp_sum<1>(n) - 2.f * bbp_sum<4>(n) - bbp_sum<5>(n) - bbp_sum<6>(n);

    // extract nth digit
    const float scaled_result = bbp_result - (int)(bbp_result) + 1.f;
    return (uint)(16.f * (scaled_result - (int)(scaled_result)));
}

static char hex_digits[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

int main()
{
    //bbp_naive(10);

    std::vector<uint> tests = {1000, 2000, 3000, 4000,5000, 10000, 50000};
    for (const auto &N : tests)
    {
        std::string decimals(N, 'x');
        auto start = std::chrono::system_clock::now();

        for (uint i = 0; i < N; ++i)
        {
            decimals[i] = hex_digits[bbp(i)];
        }
        auto end = std::chrono::system_clock::now();

    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    
    printf("%d decimals in %lld ms - %f /s\n", N, elapsed.count(), (float)1000*N /(elapsed.count()));
    }
}
