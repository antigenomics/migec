#include "migec/whitelist.hpp"

#include <zlib.h>

#include <algorithm>
#include <cmath>

#include "migec/types.hpp"

namespace migec {

Whitelist Whitelist::load(const std::string& path) {
    Whitelist w;
    gzFile fh = gzopen(path.c_str(), "rb");
    if (!fh) throw MigecError("whitelist: cannot open " + path);
    std::string line;
    char buf[1 << 16];
    size_t n_lines = 0;
    while (gzgets(fh, buf, sizeof(buf))) {
        line.assign(buf);
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
        // 10x ships its lists with a "-1" gem-group suffix on some releases.
        const size_t dash = line.find('-');
        if (dash != std::string::npos) line.resize(dash);
        if (line.empty() || line[0] == '#') continue;
        ++n_lines;
        if (w.length_ == 0) {
            w.length_ = static_cast<int>(line.size());
            if (w.length_ > kMaxBarcodeLen) {
                gzclose(fh);
                throw MigecError("whitelist: " + std::to_string(w.length_) +
                                 " nt entries do not fit the packed key (max " +
                                 std::to_string(kMaxBarcodeLen) + ")");
            }
        } else if (static_cast<int>(line.size()) != w.length_) {
            gzclose(fh);
            throw MigecError("whitelist: line " + std::to_string(n_lines) + " of " + path +
                             " is " + std::to_string(line.size()) + " nt where the file started " +
                             "with " + std::to_string(w.length_) +
                             " -- a whitelist with two lengths cannot be matched against anything");
        }
        w.keys_.push_back(pack_barcode(line));
    }
    gzclose(fh);
    if (w.keys_.empty()) throw MigecError("whitelist: " + path + " has no entries");
    std::sort(w.keys_.begin(), w.keys_.end());
    w.keys_.erase(std::unique(w.keys_.begin(), w.keys_.end()), w.keys_.end());
    return w;
}

size_t Whitelist::index_of(uint64_t key) const {
    auto it = std::lower_bound(keys_.begin(), keys_.end(), key);
    if (it == keys_.end() || *it != key) return static_cast<size_t>(-1);
    return static_cast<size_t>(it - keys_.begin());
}

bool Whitelist::contains(uint64_t key) const {
    return index_of(key) != static_cast<size_t>(-1);
}

double Whitelist::measure_background(uint64_t far_reads, uint64_t total_reads,
                                     uint64_t off_barcodes) {
    if (!total_reads || !off_barcodes) return 0.0;
    // Barcodes at distance >= 2 from every entry are off-list beyond argument, so their read
    // share is a lower bound on how much of the library is genuinely not on the list. Spread over
    // the distinct off-list barcodes, that is the prior for any one of them.
    const double share = static_cast<double>(far_reads) / static_cast<double>(total_reads);
    return share / static_cast<double>(off_barcodes);
}

std::string Whitelist::correct(std::string_view observed, std::string_view qual,
                               const std::vector<uint32_t>& counts,
                               const WhitelistParams& params, double background_prior) const {
    if (static_cast<int>(observed.size()) != length_) return {};

    // Per position: the probability the instrument got this base wrong. An N is not a failure to
    // be discarded -- it is a base consistent with all four, which is e = 0.75.
    std::vector<double> err(observed.size(), params.default_error);
    for (size_t j = 0; j < observed.size(); ++j) {
        if (base_code(observed[j]) == kInvalidBase) {
            err[j] = 0.75;
        } else if (j < qual.size()) {
            err[j] = std::clamp(phred_error(phred_from_char(qual[j])), 1e-7, 0.75);
        }
    }

    // The background: this barcode is not on the list, and every base was read correctly.
    double no_error = 1.0;
    for (double e : err) no_error *= (1.0 - e);
    double best_post = background_prior * no_error;
    double total = best_post;
    std::string best;

    // Candidates: every whitelist entry one substitution away. 3L lookups against a sorted array,
    // never a scan of the list -- 737,000 entries would be a scan per barcode.
    const uint64_t key = pack_barcode(observed);
    double sum_counts = 0.0;
    for (uint32_t c : counts) sum_counts += c;
    const double denom = sum_counts + static_cast<double>(keys_.size());

    for (int j = 0; j < length_; ++j) {
        const int shift = 62 - 2 * j;
        const uint64_t cur = (key >> shift) & 3u;
        for (uint64_t b = 0; b < 4; ++b) {
            if (b == cur) continue;
            const uint64_t cand = (key & ~(uint64_t{3} << shift)) | (b << shift);
            const size_t at = index_of(cand);
            if (at == static_cast<size_t>(-1)) continue;
            // Add-one prior on the whitelist entry's own usage: a barcode this run leans on is a
            // likelier source than one it has never seen, and the +1 keeps an unseen entry
            // possible rather than impossible.
            const double prior =
                (static_cast<double>(at < counts.size() ? counts[at] : 0u) + 1.0) / denom;
            // ...times the chance of exactly this miscall: wrong at j, right everywhere else.
            double like = err[static_cast<size_t>(j)] / 3.0;
            for (int k = 0; k < length_; ++k) {
                if (k != j) like *= (1.0 - err[static_cast<size_t>(k)]);
            }
            const double post = prior * like * (1.0 - background_prior);
            total += post;
            if (post > best_post) {
                best_post = post;
                best = unpack_barcode(cand, length_);
            }
        }
    }
    if (best.empty() || total <= 0.0) return {};
    return (best_post / total) >= params.min_posterior ? best : std::string();
}

}  // namespace migec
