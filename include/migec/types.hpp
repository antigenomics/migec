// Shared primitives: errors, base/IUPAC lookup tables, 2-bit barcode packing, Phred conversion.
//
// All lookup tables are std::array<..., 256> built once by a constexpr-ish initialiser and stored
// const, so the hot paths index them without a branch.

#ifndef MIGEC_TYPES_HPP
#define MIGEC_TYPES_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace migec {

// Every error this library raises. pybind11 maps it to RuntimeError. Messages are prefixed with
// the function name, arda-style: "mig_reader: truncated block at offset 4096".
class MigecError : public std::runtime_error {
public:
    explicit MigecError(const std::string& what) : std::runtime_error(what) {}
};

// ---------------------------------------------------------------------------------------------
// Bases

inline constexpr uint8_t kInvalidBase = 0xFF;

// ACGT -> 0..3 (upper and lower case), everything else -> kInvalidBase.
inline const std::array<uint8_t, 256>& base_code_table() {
    static const std::array<uint8_t, 256> t = [] {
        std::array<uint8_t, 256> a{};
        a.fill(kInvalidBase);
        a[static_cast<uint8_t>('A')] = a[static_cast<uint8_t>('a')] = 0;
        a[static_cast<uint8_t>('C')] = a[static_cast<uint8_t>('c')] = 1;
        a[static_cast<uint8_t>('G')] = a[static_cast<uint8_t>('g')] = 2;
        a[static_cast<uint8_t>('T')] = a[static_cast<uint8_t>('t')] = 3;
        return a;
    }();
    return t;
}

inline uint8_t base_code(char c) { return base_code_table()[static_cast<uint8_t>(c)]; }

inline char base_char(uint8_t code) { return "ACGT"[code & 3u]; }

// Complement, case-preserving; anything not ACGT maps to N.
inline const std::array<char, 256>& complement_table() {
    static const std::array<char, 256> t = [] {
        std::array<char, 256> a{};
        a.fill('N');
        const char* from = "ACGTacgtNn";
        const char* to = "TGCAtgcaNn";
        for (int i = 0; from[i]; ++i) a[static_cast<uint8_t>(from[i])] = to[i];
        return a;
    }();
    return t;
}

// Reverse-complement in place is not offered: the quality string must be reversed at the same
// time and forgetting that is the classic bug, so the two always travel together.
void reverse_complement(std::string& seq, std::string& qual);

// ---------------------------------------------------------------------------------------------
// IUPAC. Bit b of the mask is set if base code b is a member of the ambiguity set.

inline const std::array<uint8_t, 256>& iupac_mask_table() {
    static const std::array<uint8_t, 256> t = [] {
        std::array<uint8_t, 256> a{};
        auto set = [&a](char c, uint8_t m) {
            a[static_cast<uint8_t>(c)] = m;
            a[static_cast<uint8_t>(c + 32)] = m;  // lower case
        };
        set('A', 0b0001); set('C', 0b0010); set('G', 0b0100); set('T', 0b1000);
        set('U', 0b1000);
        set('R', 0b0101); set('Y', 0b1010); set('S', 0b0110); set('W', 0b1001);
        set('K', 0b1100); set('M', 0b0011);
        set('B', 0b1110); set('D', 0b1101); set('H', 0b1011); set('V', 0b0111);
        set('N', 0b1111);
        return a;
    }();
    return t;
}

inline uint8_t iupac_mask(char c) { return iupac_mask_table()[static_cast<uint8_t>(c)]; }

inline int iupac_size(uint8_t mask) {
    return ((mask >> 0) & 1) + ((mask >> 1) & 1) + ((mask >> 2) & 1) + ((mask >> 3) & 1);
}

// ---------------------------------------------------------------------------------------------
// Phred. Sanger offset 33 throughout; Illumina has not emitted offset-64 since 1.8 (2011) and
// supporting both silently is how a run comes out with every quality wrong by 31.

inline constexpr uint8_t kPhredOffset = 33;
inline constexpr uint8_t kMaxPhred = 60;

inline uint8_t phred_from_char(char c) {
    int q = static_cast<int>(static_cast<uint8_t>(c)) - kPhredOffset;
    if (q < 0) q = 0;
    if (q > kMaxPhred) q = kMaxPhred;
    return static_cast<uint8_t>(q);
}

inline char char_from_phred(uint8_t q) {
    if (q > kMaxPhred) q = kMaxPhred;
    return static_cast<char>(q + kPhredOffset);
}

// 10^(-q/10), tabulated. Callers that have a measured calibration table should use it instead;
// this is the nominal fallback.
inline const std::array<double, kMaxPhred + 1>& phred_error_table() {
    static const std::array<double, kMaxPhred + 1> t = [] {
        std::array<double, kMaxPhred + 1> a{};
        for (int q = 0; q <= kMaxPhred; ++q) a[q] = std::pow(10.0, -q / 10.0);
        return a;
    }();
    return t;
}

inline double phred_error(uint8_t q) { return phred_error_table()[q > kMaxPhred ? kMaxPhred : q]; }

// ---------------------------------------------------------------------------------------------
// Barcode packing. 2 bits per base, base 0 in the low bits, so lexicographic order of the
// unpacked string equals numeric order of the packed value only if you pack MSB-first -- we pack
// base 0 into the HIGH bits for exactly that reason: the range partition and the sort both depend
// on packed order matching barcode order.

inline constexpr int kMaxBarcodeLen = 32;  // 32 bases * 2 bits = 64

// Packs up to 32 bases. `has_n` is set if any base was not ACGT (those become A).
uint64_t pack_barcode(std::string_view seq, bool* has_n = nullptr);

std::string unpack_barcode(uint64_t packed, int len);

// Rotates a packed barcode left by `r` bases inside its `len`-base field, so the bases the range
// partition sees are a different `r` of them. This is what makes a pairwise algorithm bucketable:
// a plain range partition splits a barcode from its 1-substitution neighbour whenever the
// substitution lands in the partitioned prefix, and rotating the prefix out of the way lets the
// pair meet in a second pass. Nothing below the field moves, so an unused-tail key stays valid.
inline uint64_t rotate_barcode(uint64_t key, int len, int r) {
    if (len <= 0 || len > kMaxBarcodeLen) return key;
    r %= len;
    if (r < 0) r += len;
    if (r == 0) return key;
    const int w = 2 * len;   // significant bits, packed against the top of the word
    const int s = 64 - w;    // the unused tail below them
    const uint64_t field = key >> s;
    const uint64_t mask = w == 64 ? ~uint64_t{0} : (uint64_t{1} << w) - 1;
    return (((field << (2 * r)) | (field >> (w - 2 * r))) & mask) << s;
}

// Bucket index for the range partition: the top `bits` bits of the key. Barcodes are close to
// uniform over their alphabet, so this balances as well as a hash while preserving order.
inline uint32_t bucket_of(uint64_t key, int bits) {
    if (bits <= 0) return 0;
    return static_cast<uint32_t>(key >> (64 - bits));
}

}  // namespace migec

#endif  // MIGEC_TYPES_HPP
