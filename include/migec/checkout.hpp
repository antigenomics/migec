// checkout: find the barcode pattern in each read, extract sample/UMI, trim, and hand the barcode
// to the read header so downstream tools can see it.

#ifndef MIGEC_CHECKOUT_HPP
#define MIGEC_CHECKOUT_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "migec/pattern.hpp"

namespace migec {

enum class TrimMode {
    // Keep the read as it is. The UMI is in the header; the payload still carries the adapter.
    kNone,
    // Drop everything up to and including the matched pattern -- adapter, sample tag and UMI.
    // This is what you want before alignment: the tag is synthetic sequence and will not map.
    kPattern,
    // Drop the pattern but keep whatever preceded it. Rarely what you want; provided because some
    // chemistries put real sequence 5' of the tag.
    kPatternOnly,
};

struct CheckoutParams {
    MatchParams match;
    TrimMode trim = TrimMode::kPattern;
    // Reads shorter than this after trimming are dropped as uninformative.
    int min_payload = 1;
    // Also require the UMI to have no ambiguous base. Off by default: an N in a UMI is a reason
    // to be less certain, not a reason to throw the molecule away.
    bool reject_umi_n = false;
    // Minimum Phred across the UMI bases. 0 disables. MIGEC used 15 and MAGERI 20, both as hard
    // drops; the default here is 0 because discarding a read for a low-quality UMI base loses
    // sequence information that the correction step can often recover.
    int min_umi_quality = 0;
};

struct CheckoutCounters {
    uint64_t total = 0;
    uint64_t assigned = 0;
    uint64_t unmatched = 0;
    uint64_t ambiguous = 0;
    uint64_t short_payload = 0;
    uint64_t bad_umi = 0;
    std::vector<uint64_t> per_sample;
};

struct CheckoutRead {
    bool ok = false;
    int sample = -1;
    std::string umi;
    std::string umi_qual;
    // Views into the input read, after trimming.
    std::string_view seq;
    std::string_view qual;
    double score = 0.0;
};

// Stateless apart from the pattern set and the counters.
class Checkout {
public:
    Checkout(const PatternSet& patterns, CheckoutParams params);

    CheckoutRead process(std::string_view seq, std::string_view qual);

    const CheckoutCounters& counters() const { return counters_; }

    // The SAM-style comment to append to a FASTQ header. TAB-separated, because `bwa mem -C` and
    // `minimap2 -y` copy the comment verbatim into the SAM record and it has to be conformant
    // there. `sample` may be empty.
    static std::string header_tags(const std::string& umi, const std::string& umi_qual,
                                   const std::string& sample);

private:
    const PatternSet& patterns_;
    CheckoutParams params_;
    CheckoutCounters counters_;
};

}  // namespace migec

#endif  // MIGEC_CHECKOUT_HPP
