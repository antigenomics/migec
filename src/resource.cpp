#include "migec/resource.hpp"

#include <sys/resource.h>

#include <thread>

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
    const unsigned n = std::thread::hardware_concurrency();
    return n ? n : 1u;
}

}  // namespace migec
