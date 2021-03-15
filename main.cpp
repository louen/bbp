#include "CoreCpp/CoreMacros.hpp"
#include "CoreCpp/CoreStrings.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

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
    static const char *ansi_white = "\033[0;37m";
    static const char *ansi_bold_green = "\033[1;32m";

    printf(ansi_white);

    int x, y, n;
    unsigned char *data = stbi_load("pi.png", &x, &y, &n, 0);

    uint index = 0;

    std::vector<uint> digits(x * y);

    for (int i = 20; i < x - 20; ++i)
    {
        for (int j = 0; j < y; ++j)
        {
            if (index == 1)
            {
                printf(".");
            }
            else
            {
                printf("%c", hex_digits[bbp(index ? index - 1 : index)]);
            }
            ++index;

            const uint pixel = (i * x + j) * n;
            if (j < y - 1)
            {
                const uint next_pixel = (i * x + j + 1) * n;

                if (data[pixel] == 0 && data[next_pixel] > 0)
                {
                    printf(ansi_bold_green);
                }
                else if (data[pixel] > 0 && data[next_pixel] == 0)
                {
                    printf(ansi_white);
                }
            }
            else
            {
                printf(ansi_white);
            }
        }
        printf("\n");
    }
}
