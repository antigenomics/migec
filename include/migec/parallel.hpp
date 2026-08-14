// One thread helper, shared by every stage that has work to spread.
//
// Never: nothing may throw out of a worker thread. An escaping exception is std::terminate --
// SIGABRT, no message, no flush, no output. Workers capture, the caller rethrows, and the first
// exception by INDEX wins rather than the first by wall clock, so a failure reports the same way
// whatever the thread count.
//
// Never: the work is indexed, never queued, and results go into per-index slots the caller owns.
// A worker never appends to a shared container and never takes a lock on the hot path, which is
// what makes "the output does not depend on -t" a property of the shape rather than a promise.

#ifndef MIGEC_PARALLEL_HPP
#define MIGEC_PARALLEL_HPP

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

#include "migec/resource.hpp"

namespace migec {

// Threads to use for `items` units of work: `requested`, or one per core when it is 0, never more
// than there is work for and never fewer than one.
inline int worker_count(int requested, size_t items) {
    int n = requested > 0 ? requested : static_cast<int>(hardware_threads());
    if (n < 1) n = 1;
    if (items && static_cast<size_t>(n) > items) n = static_cast<int>(items);
    return n;
}

// Runs `fn(i, worker)` for every i in [0, items), on `threads` threads. `worker` is a stable index
// in [0, threads) so a caller can give each thread its own scratch without a lock.
//
// Work is claimed with an atomic counter rather than split into contiguous blocks: the per-item
// cost here is wildly uneven -- one barcode's neighbourhood is empty and the next one's is full,
// one bucket holds ten molecules and the next ten million -- and a static split would leave every
// thread but one idle at the end.
//
// Never: the counter is claimed in BATCHES, not one item at a time. Measured on assemble's
// partition, where an item is one read's tag scan: 21% of all CPU samples across every thread sat
// on that single `ldadd` instruction, more than the parse it was handing out. Sixteen cores
// hammering one cache line is the work, not the work. The batch is sized so that each worker takes
// ~8 turns -- enough turns for the uneven case to even out, few enough atomics to disappear -- and
// collapses to 1 when there are few items, which is exactly the uneven case (one bucket per item)
// that the counter exists for.
template <typename Fn>
void parallel_for(size_t items, int threads, Fn&& fn) {
    const int n = worker_count(threads, items);
    if (items == 0) return;
    if (n == 1) {
        for (size_t i = 0; i < items; ++i) fn(i, 0);
        return;
    }

    // Capped, so a ten-million-item scan still hands out ten thousand batches rather than eight:
    // the tail imbalance of a batch is paid by whichever worker draws the last one.
    const size_t grab =
        std::min<size_t>(1024, std::max<size_t>(1, items / (static_cast<size_t>(n) * 8)));
    std::atomic<size_t> next{0};
    std::mutex err_mutex;
    std::exception_ptr err;
    size_t err_at = items;

    auto run = [&](int worker) {
        for (;;) {
            const size_t first = next.fetch_add(grab, std::memory_order_relaxed);
            if (first >= items) return;
            const size_t last = std::min(items, first + grab);
            for (size_t i = first; i < last; ++i) {
                try {
                    fn(i, worker);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(err_mutex);
                    if (i < err_at) {  // lowest index wins: the message does not depend on timing
                        err_at = i;
                        err = std::current_exception();
                    }
                    return;
                }
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(n) - 1);
    // Never: a failed SPAWN must not abort. `emplace_back` throws std::system_error under a
    // thread or RLIMIT_NPROC cap, and unwinding then runs ~thread() on the threads already
    // started -- which is std::terminate: SIGABRT, no message, nothing flushed, which is the one
    // outcome this file exists to prevent. The threads that did start finish their work, and the
    // rest of the items are run right here, so the answer is the same and only the wall clock
    // changes.
    std::exception_ptr spawn_err;
    for (int t = 1; t < n; ++t) {
        try {
            pool.emplace_back([&run, t] { run(t); });
        } catch (...) {
            spawn_err = std::current_exception();
            break;
        }
    }
    run(0);
    for (std::thread& th : pool) th.join();
    if (err) std::rethrow_exception(err);
    // The work is done -- `run` claims batches until the queue is empty, so whoever ran took the
    // rest -- and the spawn failure is reported only if it left something undone.
    if (spawn_err && next.load(std::memory_order_relaxed) < items) {
        std::rethrow_exception(spawn_err);
    }
}

}  // namespace migec

#endif  // MIGEC_PARALLEL_HPP
