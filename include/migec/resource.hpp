// Wall clock and peak memory. Both are reported on every run rather than hidden behind a verbose
// flag: throughput and resident size are the two numbers that decide whether a pipeline can be run
// at all, and a tool that only prints them when asked gets run once without asking.

#ifndef MIGEC_RESOURCE_HPP
#define MIGEC_RESOURCE_HPP

#include <chrono>
#include <cstddef>

namespace migec {

// Peak resident set size of this process in bytes. 0 when the platform will not say.
size_t peak_rss_bytes();

// Number of hardware threads, never 0.
unsigned hardware_threads();

class Stopwatch {
public:
    Stopwatch() : start_(std::chrono::steady_clock::now()) {}
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }

private:
    std::chrono::steady_clock::time_point start_;
};

}  // namespace migec

#endif  // MIGEC_RESOURCE_HPP
