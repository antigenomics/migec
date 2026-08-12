#include "migec/types.hpp"

#include <algorithm>

namespace migec {

void reverse_complement(std::string& seq, std::string& qual) {
    if (!qual.empty() && qual.size() != seq.size()) {
        throw MigecError("reverse_complement: sequence and quality lengths differ (" +
                         std::to_string(seq.size()) + " vs " + std::to_string(qual.size()) + ")");
    }
    const auto& comp = complement_table();
    const size_t n = seq.size();
    for (size_t i = 0, j = n ? n - 1 : 0; i < j; ++i, --j) {
        char a = comp[static_cast<uint8_t>(seq[i])];
        char b = comp[static_cast<uint8_t>(seq[j])];
        seq[i] = b;
        seq[j] = a;
    }
    if (n & 1u) seq[n / 2] = comp[static_cast<uint8_t>(seq[n / 2])];
    std::reverse(qual.begin(), qual.end());
}

uint64_t pack_barcode(std::string_view seq, bool* has_n) {
    if (seq.size() > static_cast<size_t>(kMaxBarcodeLen)) {
        throw MigecError("pack_barcode: barcode longer than 32 bases (" +
                         std::to_string(seq.size()) + ")");
    }
    uint64_t packed = 0;
    bool n = false;
    // Base 0 goes into the high bits so that packed order == lexicographic order of the barcode,
    // which is what makes the range partition and the sort agree.
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t c = base_code(seq[i]);
        if (c == kInvalidBase) {
            n = true;
            c = 0;  // an N is stored as A; the flag carries the ambiguity
        }
        packed |= static_cast<uint64_t>(c) << (62 - 2 * i);
    }
    if (has_n) *has_n = n;
    return packed;
}

std::string unpack_barcode(uint64_t packed, int len) {
    if (len < 0 || len > kMaxBarcodeLen) {
        throw MigecError("unpack_barcode: length out of range (" + std::to_string(len) + ")");
    }
    std::string out(static_cast<size_t>(len), 'A');
    for (int i = 0; i < len; ++i) {
        out[static_cast<size_t>(i)] = base_char(static_cast<uint8_t>((packed >> (62 - 2 * i)) & 3u));
    }
    return out;
}

}  // namespace migec
