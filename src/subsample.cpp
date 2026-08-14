#include "migec/subsample.hpp"

#include <algorithm>
#include <unordered_map>

#include "migec/fastq.hpp"
#include "migec/resource.hpp"
#include "migec/types.hpp"

namespace migec {

SubsampleStats subsample(const SubsampleRequest& request) {
    Stopwatch clock;
    SubsampleStats stats;
    if (request.per_10k == 0 || request.per_10k > 10000) {
        throw MigecError("subsample: keep must be in 1..10000 ten-thousandths");
    }

    FastqWriter writer(request.output, request.gzip_level);
    FastqReader reader(request.input);
    FastqRecord rec;
    // Only over the barcodes that were KEPT, so this is bounded by the fixture and not by the
    // library it came from -- which is the whole point of running this before anything else.
    // Counted, not just seen: the median and the examples both need the per-barcode depth, and
    // a mean alone cannot tell a fixture that kept the distribution from one that flattened it.
    std::unordered_map<uint64_t, uint64_t> kept_keys;
    size_t key_length = 0;
    while (reader.next(rec)) {
        ++stats.reads;
        const std::string_view umi = tag_value(rec.comment, "RX:Z:");
        if (umi.empty()) { ++stats.reads_without_umi; continue; }
        // Select on the CELL when there is one: a kept cell keeps every molecule in it. Hashing
        // cell+umi together would sample molecules independently, which is the read-sampling
        // mistake wearing a different hat.
        const std::string_view cell = request.by_cell ? tag_value(rec.comment, "CB:Z:")
                                                      : std::string_view{};
        const uint64_t selector = cell.empty() ? pack_barcode(umi) : pack_barcode(cell);
        if (!keeps(selector, request.per_10k)) continue;
        ++stats.reads_kept;
        // Counted as MOLECULES, whatever the selection was on, because reads-per-barcode is the
        // number that says whether the size distribution survived.
        if (cell.empty()) {
            ++kept_keys[selector];
            key_length = umi.size();
        } else {
            std::string joined;
            joined.reserve(cell.size() + umi.size());
            joined.append(cell);
            joined.append(umi);
            ++kept_keys[pack_barcode(joined)];
            key_length = joined.size();
        }
        writer.write(rec.name, rec.comment, rec.seq, rec.qual);
    }
    writer.close();
    stats.barcodes_seen = kept_keys.size();

    if (!kept_keys.empty()) {
        std::vector<uint64_t> keys;
        std::vector<uint64_t> depths;
        keys.reserve(kept_keys.size());
        depths.reserve(kept_keys.size());
        for (const auto& [key, count] : kept_keys) {
            keys.push_back(key);
            depths.push_back(count);
        }
        const size_t mid = depths.size() / 2;
        std::nth_element(depths.begin(), depths.begin() + static_cast<ptrdiff_t>(mid),
                         depths.end());
        stats.reads_per_barcode_median = depths[mid];
        stats.reads_per_barcode_max = *std::max_element(depths.begin(), depths.end());
        // Key order, so which barcodes are shown does not depend on how deeply they were
        // sequenced -- see the note on SubsampleStats::examples.
        constexpr size_t kExamples = 5;
        const size_t shown = std::min(kExamples, keys.size());
        std::partial_sort(keys.begin(), keys.begin() + static_cast<ptrdiff_t>(shown), keys.end());
        for (size_t i = 0; i < shown; ++i) {
            stats.examples.emplace_back(unpack_barcode(keys[i], static_cast<int>(key_length)),
                                        kept_keys[keys[i]]);
        }
    }
    stats.wall_seconds = clock.seconds();
    stats.peak_rss_bytes = peak_rss_bytes();
    return stats;
}

}  // namespace migec
