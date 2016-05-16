#include "amplicon_index.h"

namespace atlas {

amplicon_index::amplicon_index(amplicon_pairs const& amplicons, size_t k) {
    k_ = k;
    modulo_ = 1;
    for (size_t i = 1; i < k_; ++i) {
        modulo_ *= 4;
    }

    int amplicon_number = 1;
    for (auto const& amplicon : amplicons) {
        add_amplicon(amplicon.first, amplicon_number);
        add_amplicon(amplicon.second, -amplicon_number);
        ++amplicon_number;
    }
}

amplicon_index::amplicon_index(std::istream& in) {
    in >> k_;
    modulo_ = 1;
    for (size_t i = 1; i < k_; ++i) {
        modulo_ *= 4;
    }

    while (true) {
        int kmer;
        if (!(in >> kmer)) {
            return;
        }

        std::vector<int> amplicons;
        int amplicon_number;
        in >> amplicon_number;
        while (amplicon_number) {
            amplicons.push_back(amplicon_number);
            in >> amplicon_number;
        }

        dict_[kmer] = std::move(amplicons);
    }
}

void amplicon_index::add_amplicon(std::string const& seq, int amplicon_number) {
    u64 kmer = 0;

    for (size_t i = 0; i < k_ - 1; ++i) {
        kmer = 4 * kmer + get_nt_number(seq[i]);
    }

    for (size_t i = k_ - 1; i < seq.size(); ++i) {
        kmer = 4 * (kmer % modulo_) + get_nt_number(seq[i]);
        dict_[kmer].push_back(amplicon_number);
    }
}

bool amplicon_index::contains(u64 kmer) const {
    return dict_.find(kmer) != dict_.end();
}

std::vector<int> const& amplicon_index::get_all_amplicons(u64 kmer) const {
    return dict_.find(kmer)->second;
}

size_t amplicon_index::get_k() const {
    return k_;
}

void amplicon_index::save(std::ostream& out) const {
    out << k_ << "\n";
    for (auto const& entry : dict_) {
        out << entry.first;
        for (int amplicon : entry.second) {
            out << " " << amplicon;
        }
        out << " 0" << std::endl;
    }
}

amplicon_pairs create_amplicon_pairs(std::istream& amplicons_stream) {
    amplicon_pairs pairs;
    for (read const& r : read_fasta(amplicons_stream)) {
        pairs.push_back(std::make_pair(r.get_seq(), reverse_complement(r.get_seq())));
    }
    return pairs;
}

} // namespace atlas
