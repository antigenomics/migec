#pragma once

#include "migec/parallel.hpp"
#include "migec/types.hpp"
#include "migec/umi_stats.hpp"
#include "seqtree/seqtree.hpp"

namespace migec::detail {

// The same native batch API owns distance-one search and broader edit candidates.
// Reference ids stay aligned with the caller-owned barcode/evidence tables.
inline std::vector<std::vector<seqtree::Hit>> barcode_neighbors(
    const std::vector<UmiCounts::Entry>& entries, int length, int substitutions, int threads) {
    std::vector<std::string> sequences;
    sequences.reserve(entries.size());
    for (const auto& e : entries) sequences.push_back(unpack_barcode(e.key, length));
    auto index = seqtree::Index::build(sequences, seqtree::Alphabet::Nucleotide);
    seqtree::SearchParams scope;
    scope.engine = seqtree::Engine::SeqTm;
    scope.max_substitutions = static_cast<uint16_t>(substitutions);
    return index->search_batch(sequences, scope, worker_count(threads, entries.size()));
}

}  // namespace migec::detail
