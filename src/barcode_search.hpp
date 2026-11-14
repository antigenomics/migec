#pragma once

#include <limits>

#include "migec/parallel.hpp"
#include "migec/types.hpp"
#include "migec/umi_stats.hpp"
#include "seqtree/seqtree.hpp"

namespace migec::detail {

// One owned native index per caller table; no persisted results or second search engine.
struct BarcodeSearch {
    std::vector<std::string> queries;
    std::unique_ptr<seqtree::Index> index;

    BarcodeSearch(const std::vector<UmiCounts::Entry>& entries, int length) {
        if (length < 0 || length > 32) throw MigecError("invalid barcode search length");
        queries.reserve(entries.size());
        uint64_t previous = 0;
        bool first = true;
        for (const auto& entry : entries) {
            if ((length == 0 && entry.key != 0) || (!first && entry.key <= previous) ||
                (length > 0 && length < 32 && (entry.key << (2 * length)) != 0))
                throw MigecError("barcode search requires sorted unique keys with zero padding");
            first = false;
            previous = entry.key;
            queries.push_back(unpack_barcode(entry.key, length));
        }
        index = seqtree::Index::build(queries, seqtree::Alphabet::Nucleotide);
    }

    void consume(int substitutions, int threads,
                 const std::function<void(size_t, size_t, seqtree::CountBuffer&&)>& callback) const {
        if (substitutions < 0 || substitutions > 2)
            throw MigecError("barcode search consumer supports at most two substitutions");
        seqtree::SearchParams scope;
        scope.engine = seqtree::Engine::SeqTm;
        scope.max_substitutions = static_cast<uint16_t>(substitutions);
        constexpr size_t chunk = 32;
        // Unique DNA tags of at most 32 nt have at most 1+3*32+9*C(32,2)=4561
        // radius-two hits. A 32-query primitive block therefore occupies <4 MiB,
        // independent of table size. The generic API's conservative index-size
        // ceiling is supplied separately so it does not shrink sparse batches to one query.
        if (queries.size() > std::numeric_limits<size_t>::max() / (chunk * 6 * sizeof(uint32_t)))
            throw MigecError("barcode search buffer ceiling exceeds the address space");
        const size_t ceiling = std::max<size_t>(1, queries.size()) * chunk * 6 * sizeof(uint32_t);
        index->search_batch_consume(queries, scope, callback,
                                   worker_count(threads, queries.size()), chunk, ceiling);
    }
};

inline std::vector<std::vector<seqtree::Hit>> barcode_neighbors(
    const std::vector<UmiCounts::Entry>& entries, int length, int substitutions, int threads) {
    BarcodeSearch search(entries, length);
    seqtree::SearchParams scope;
    scope.engine = seqtree::Engine::SeqTm;
    scope.max_substitutions = static_cast<uint16_t>(substitutions);
    return search.index->search_batch(search.queries, scope, worker_count(threads, entries.size()));
}

}  // namespace migec::detail
