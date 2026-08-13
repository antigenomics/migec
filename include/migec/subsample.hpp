// subsample: take ALL the reads of a fraction of the UMIs.
//
// ⛔ Never take a fraction of the READS. Ten thousand random reads of a library sequenced at four
// reads per molecule gives ten thousand molecules seen once each, which is not a smaller version
// of the library -- it is a different library, with the MIG size distribution destroyed and every
// consensus reduced to a single read. Every fixture built that way silently tests nothing.
//
// ⛔ Nor the FIRST N distinct UMIs. A UMI with 100 reads appears in the first thousand reads about
// a hundred times more often than a singleton does, so first-appearance order oversamples large
// MIGs -- again destroying the distribution the fixture exists to show.
//
// So: hash the barcode, keep the read if the hash falls in the kept range. One streaming pass, no
// sort, no memory, and the same reads chosen every time on any machine.

#ifndef MIGEC_SUBSAMPLE_HPP
#define MIGEC_SUBSAMPLE_HPP

#include <cstdint>
#include <string>
#include <string_view>

namespace migec {

// splitmix64, on the packed barcode. Not blake2b: this needs no dependency, is one multiply-xor
// chain rather than a compression function, and the only property being asked of it is that the
// low bits are uncorrelated with barcode content. It is written out here rather than called from
// a library precisely so that the selection can be reproduced by anything, in any language.
inline uint64_t barcode_hash(uint64_t key) {
    uint64_t z = key + 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// The barcode is kept when `barcode_hash(key) % 10000 < per_10k`.
inline bool keeps(uint64_t key, uint32_t per_10k) {
    return barcode_hash(key) % 10000ull < per_10k;
}

struct SubsampleRequest {
    std::string input;
    std::string output;
    // Barcodes kept, in ten-thousandths. 100 = 1%.
    uint32_t per_10k = 100;
    // Select on the CELL barcode when the reads carry one, so a kept cell keeps every molecule in
    // it. Hashing cell+umi together would sample molecules independently and give a fixture of
    // thousands of cells holding a handful of molecules each -- the read-sampling mistake wearing
    // a different hat. Falls back to the UMI when there is no cell barcode.
    bool by_cell = true;
    int gzip_level = 6;
};

struct SubsampleStats {
    uint64_t reads = 0;
    uint64_t reads_kept = 0;
    uint64_t barcodes_seen = 0;   // distinct, among kept reads
    uint64_t reads_without_umi = 0;
    double wall_seconds = 0.0;
};

SubsampleStats subsample(const SubsampleRequest& request);

}  // namespace migec

#endif  // MIGEC_SUBSAMPLE_HPP
