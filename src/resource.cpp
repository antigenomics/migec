#include "migec/resource.hpp"

#include <sys/resource.h>

#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <fstream>
#include <string_view>
#include <thread>

#ifdef __linux__
#include <sched.h>
#endif

namespace migec {

size_t peak_rss_bytes() {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
#ifdef __APPLE__
    // Darwin reports ru_maxrss in bytes, Linux in kilobytes. Getting this wrong is a factor of
    // 1024 in the one number a user checks before deciding a run will not fit.
    return static_cast<size_t>(ru.ru_maxrss);
#else
    return static_cast<size_t>(ru.ru_maxrss) * 1024u;
#endif
}

unsigned hardware_threads() {
    unsigned count = std::max(1u, std::thread::hardware_concurrency());
#ifdef __linux__
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0)
        count = std::min(count, static_cast<unsigned>(std::max(1, CPU_COUNT(&affinity))));
    long long quota = 0, period = 0;
    std::ifstream limits("/sys/fs/cgroup/cpu.max");
    if (!(limits >> quota >> period)) {
        std::ifstream legacy_quota("/sys/fs/cgroup/cpu/cpu.cfs_quota_us");
        std::ifstream legacy_period("/sys/fs/cgroup/cpu/cpu.cfs_period_us");
        legacy_quota >> quota;
        legacy_period >> period;
    }
    if (quota > 0 && period > 0)
        count = static_cast<unsigned>(std::min<long long>(count, std::max(1LL, quota / period)));
#endif
    for (const char* name : {"SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"}) {
        if (const char* value = std::getenv(name)) {
            const std::string_view text(value);
            unsigned allocated = 0;
            const auto parsed = std::from_chars(text.data(), text.data() + text.size(), allocated);
            if (parsed.ec == std::errc{} && parsed.ptr == text.data() + text.size() && allocated)
                count = std::min(count, allocated);
        }
    }
    return count;
}

}  // namespace migec
