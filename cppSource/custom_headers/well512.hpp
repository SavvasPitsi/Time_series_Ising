#include <iostream>
#include <array>
#include <cmath>  // For std::nextafter
#include <stdint.h>

class WELL512 {
private:
    static constexpr uint32_t STATE_SIZE = 16;  // State size
    std::array<uint32_t, STATE_SIZE> state;
    uint32_t index = 0;

public:
    // Constructor: Initializes state using a seed
    WELL512(uint32_t seed = 5489U) {
        state[0] = seed;
        for (uint32_t i = 1; i < STATE_SIZE; ++i) {
            state[i] = (1812433253U * (state[i - 1] ^ (state[i - 1] >> 30)) + i);
        }
    }

    // Generate next random uint32_t
    uint32_t next() {
        uint32_t a, b, c, d;
        a = state[index];
        c = state[(index + 13) & 15];
        b = a ^ (a << 16) ^ c ^ (c << 15);
        c = state[(index + 9) & 15];
        c ^= (c >> 11);
        a = state[index] = b ^ c;
        d = a ^ ((a << 5) & 0xDA442D24U);
        index = (index + 15) & 15;
        return state[index] = state[index] ^ b ^ d ^ (state[(index + 7) & 15] ^ (state[(index + 7) & 15] >> 2));
    }

    // Generate next random double in [0, 1)
    double nextDouble() {
        return (next() >> 8) * (1.0 / 16777216.0);  // Scales to [0,1)
    }
};

// Global WELL512 generator instance
static WELL512 rng;

// Function to generate a random number in [0, max)
double random(double max) {
    return rng.nextDouble() * std::nextafter(max, 0.0);
}

