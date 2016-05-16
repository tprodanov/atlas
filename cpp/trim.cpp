#include "trim.h"

namespace atlas {

float trim_parameters::score_threshold_ratio = 0.5;
float trim_parameters::gap_penalty = -1.4;
float trim_parameters::mismatch_penalty = -1.4;
size_t trim_parameters::degenerate_nucleotides_count = 2;
bool trim_parameters::paired_local = false;

float local_align(std::string const& seq, std::string const& amplicon, size_t& l, size_t& r) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    size_t *starts = new size_t[(n + 1) * (m + 1)];
    float score = -1e6;

    for (size_t j = 0; j <= n; ++j) {
        h[j] = 0;
        starts[j] = j;
    }

    for (size_t i = 1; i <= m; ++i) {
        h[i * (n + 1)] = 0;
        starts[i * (n + 1)] = 0;

        for (size_t j = 1; j <= n; ++j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j - 1] + trim_parameters::gap_penalty;
            starts[i * (n + 1) + j] = starts[i * (n + 1) + j - 1];

            if (h[(i - 1) * (n + 1) + j] + trim_parameters::gap_penalty > h[i * (n + 1) + j]) {
                h[i * (n + 1) + j] = h[(i - 1) * (n + 1) + j] + trim_parameters::gap_penalty;
                starts[i * (n + 1) + j] = starts[(i - 1) * (n + 1) + j];
            }

            if (h[(i - 1) * (n + 1) + j - 1]
                    + (seq[j - 1] == amplicon[i - 1] ? 1.0f : trim_parameters::mismatch_penalty)
                    > h[i * (n + 1) + j]) {
                h[i * (n + 1) + j] = h[(i - 1) * (n + 1) + j - 1]
                        + (seq[j - 1] == amplicon[i - 1] ? 1.0f : trim_parameters::mismatch_penalty);
                starts[i * (n + 1) + j] = starts[(i - 1) * (n + 1) + j - 1];
            }

            if (h[i * (n + 1) + j] < 0) {
                h[i * (n + 1) + j] = 0;
                starts[i * (n + 1) + j] = j;
            }

            if (h[i * (n + 1) + j] > score) {
                score = h[i * (n + 1) + j];
                l = starts[i * (n + 1) + j];
                r = j;
            }
        }
    }

    delete[] starts;
    delete[] h;
    return score;
}

float align_from_the_left(std::string const& seq, std::string const& amplicon, size_t&, size_t& r) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    float score = -1e6;

    for (size_t i = 0; i <= m; ++i) {
        h[i * (n + 1)] = 0;

        for (size_t j = 1; j <= n; ++j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j - 1] + trim_parameters::gap_penalty;

            if (i) {
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j],
                                              h[(i - 1) * (n + 1) + j] + trim_parameters::gap_penalty);
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j],
                        h[(i - 1) * (n + 1) + j - 1] +
                            (seq[j - 1] == amplicon[i - 1] ? 1.0f : trim_parameters::mismatch_penalty));
            }

            if (h[i * (n + 1) + j] >= score) {
                score = h[i * (n + 1) + j];
                r = j;
            }
        }
    }

    delete[] h;
    return score;
}

float align_from_the_right(std::string const& seq, std::string const& amplicon, size_t& l, size_t&) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    float score = -1e6;

    for (size_t i = m; i <= m; --i) {
        h[i * (n + 1) + n] = 0;

        for (size_t j = n - 1; j <= n; --j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j + 1] + trim_parameters::gap_penalty;

            if (i < m) {
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j],
                                              h[(i + 1) * (n + 1) + j] + trim_parameters::gap_penalty);
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j],
                        h[(i + 1) * (n + 1) + j + 1] +
                            (seq[j] == amplicon[i] ? 1.0f : trim_parameters::mismatch_penalty));
            }

            if (h[i * (n + 1) + j] >= score) {
                score = h[i * (n + 1) + j];
                l = j;
            }
        }
    }

    delete[] h;
    return score;
}

std::set<int> possible_amplicons(std::string const& seq, amplicon_index const& ai) {
    size_t n = seq.size();
    size_t k = ai.get_k();
    if (n < k) {
        return std::set<int>();
    }

    static size_t const kmers_count = 7;
    size_t kmer_pos[kmers_count] = {n / 2 - k / 2,
                                    k, n - 2 * k,
                                    n / 3 - k / 2, 2 * n / 3 - k / 2,
                                    0, n - k};

    std::set<int> result;
    for (size_t i = 0; i < kmers_count; ++i) {
        if (i >= 3 && !result.empty())
            return result;
        if (kmer_pos[i] > n || kmer_pos[i] + k > n)
            continue;

        for (u64 kmer : seq_kmer(seq, k, kmer_pos[i], trim_parameters::degenerate_nucleotides_count)) {
            if (!ai.contains(kmer))
                continue;
            for (int amplicon : ai.get_all_amplicons(kmer))
                result.insert(amplicon);
        }
    }
    return result;
}

std::set<int> possible_amplicons(std::string const& seq1,
                                 std::string const& seq2,
                                 amplicon_index const& ai,
                                 bool rc) {
    int ix_multiplier = rc ? -1 : 1;

    size_t n1 = seq1.size();
    size_t n2 = seq2.size();

    size_t k = ai.get_k();
    if (n1 < k || n2 < k) {
        return std::set<int>();
    }

    static size_t const kmers_count = 7;
    size_t kmer_pos1[kmers_count] = {n1 / 2 - k / 2,
                                     k, n1 - 2 * k,
                                     n1 / 3 - k / 2, 2 * n1 / 3 - k / 2,
                                     0, n1 - k};
    size_t kmer_pos2[kmers_count] = {n2 / 2 - k / 2,
                                     k, n2 - 2 * k,
                                     n2 / 3 - k / 2, 2 * n2 / 3 - k / 2,
                                     0, n2 - k};

    std::vector<int> amplicons1;
    std::set<int> amplicons2;
    std::set<int> result;

    for (size_t i = 0; i < kmers_count; ++i) {
        if (i >= 3 && i % 2 == 1) {
            for (int ix : amplicons1)
                if (amplicons2.find(ix_multiplier * ix) != amplicons2.end()) {
                    result.insert(ix);
                }

            if (!result.empty())
                return result;
        }

        if (kmer_pos1[i] <= n1 && kmer_pos1[i] + k <= n1) {
            for (u64 kmer : seq_kmer(seq1, k, kmer_pos1[i], trim_parameters::degenerate_nucleotides_count)) {
                if (ai.contains(kmer)) {
                    for (int amplicon : ai.get_all_amplicons(kmer))
                        amplicons1.push_back(amplicon);
                }
            }
        }

        if (kmer_pos2[i] <= n2 && kmer_pos2[i] + k <= n2) {
            for (u64 kmer : seq_kmer(seq2, k, kmer_pos2[i], trim_parameters::degenerate_nucleotides_count)) {
                if (ai.contains(kmer)) {
                    for (int amplicon : ai.get_all_amplicons(kmer))
                        amplicons2.insert(amplicon);
                }
            }
        }
    }

    for (int ix : amplicons1)
        if (amplicons2.find(ix_multiplier * ix) != amplicons2.end())
            result.insert(ix);

    if (result.empty()) {
        for (int amplicon : amplicons1)
            amplicons2.insert(amplicon);
        return amplicons2;
    }

    return result;
}

bool trim_read(read& r, amplicon_index const& ai, amplicon_pairs const& amplicons) {
    std::string const& seq = r.get_seq();
    float best_score = -1;
    size_t best_begin = 0;
    size_t best_end = seq.size();

    for (int ix : possible_amplicons(seq, ai)) {
        size_t left;
        size_t right;
        float score = ix > 0
            ? local_align(seq, amplicons[ix - 1].first, left, right)
            : local_align(seq, amplicons[-ix - 1].second, left, right);
        if (score > best_score) {
            best_score = score;
            best_begin = left;
            best_end = right;
        }
    }

    if (best_score >= trim_parameters::score_threshold_ratio * seq.size()) {
        r.trim(best_begin, best_end);
        return true;
    }
    return false;
}

// rc - r2 is reverse complement comparing with r1
bool trim_reads(read& r1, read& r2, amplicon_index const& ai, amplicon_pairs const& amplicons, bool rc) {
    std::string const& seq1 = r1.get_seq();
    std::string const& seq2 = r2.get_seq();
    float best_score = -1;

    size_t begin1 = 0;
    size_t end1 = seq1.size();
    size_t begin2 = 0;
    size_t end2 = seq2.size();

    auto& align1 = trim_parameters::paired_local ? local_align : align_from_the_right;
    auto& align2 = trim_parameters::paired_local ? local_align : align_from_the_left;

    for (int ix : possible_amplicons(seq1, seq2, ai, rc)) {
        size_t left1 = 0;
        size_t right1 = seq1.size();
        float score1 = ix > 0
            ? align1(seq1, amplicons[ix - 1].first, left1, right1)
            : align1(seq1, amplicons[-ix - 1].second, left1, right1);

        size_t left2 = 0;
        size_t right2 = seq2.size();
        float score2;
        if (rc) {
            score2 = ix > 0
                ? align2(seq2, amplicons[ix - 1].second, left2, right2)
                : align2(seq2, amplicons[-ix - 1].first, left2, right2);
        } else {
            score2 = ix > 0
                ? align2(seq2, amplicons[ix - 1].first, left2, right2)
                : align2(seq2, amplicons[-ix - 1].second, left2, right2);
        }

        if (score1 + score2 > best_score) {
            best_score = score1 + score2;
            begin1 = left1;
            end1 = right1;

            begin2 = left2;
            end2 = right2;
        }
    }

    if (best_score >= trim_parameters::score_threshold_ratio * (seq1.size() + seq2.size())) {
        r1.trim(begin1, end1, r2, begin2, end2);
        return true;
    }
    return false;
}

void trim_and_write(read& r, amplicon_index const& ai, amplicon_pairs const& amplicons,
                    std::ostream& out, std::ostream& err) {
    try {
        if (trim_read(r, ai, amplicons)) {
            r.write_to(out);
        } else {
            r.write_to(err);
        }
    } catch (unexpected_nt const&) {
        r.write_to(err);
    }
}

// rc - r2 is reverse complement comparing with r1
void trim_and_write(read& r1, read& r2, bool rc, amplicon_index const& ai, amplicon_pairs const& amplicons,
                    std::ostream& out1, std::ostream& out2, std::ostream& err1, std::ostream& err2) {
    try {
        if (trim_reads(r1, r2, ai, amplicons, rc)) {
            r1.write_to(out1);
            r2.write_to(out2);
        } else {
            r1.write_to(err1);
            r2.write_to(err2);
        }
    } catch (unexpected_nt const&) {
        r1.write_to(err1);
        r2.write_to(err2);
    }
}

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq, std::ostream& output, std::ostream& error) {
    fastq_read r;
    while (fastq_read::from_stream(fastq, r)) {
        trim_and_write(r, ai, amplicons, output, error);
    }
}

void trim_fastq(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq1,  std::istream& fastq2,
          std::ostream& output1, std::ostream& output2,
          std::ostream& error1,  std::ostream& error2) {
    fastq_read r1;
    fastq_read r2;

    while (fastq_read::from_stream(fastq1, r1) && fastq_read::from_stream(fastq2, r2)) {
        trim_and_write(r1, r2, true, ai, amplicons, output1, output2, error1, error2);
    }
}

void trim_sam(amplicon_index const& ai, amplicon_pairs const& amplicons,
              std::istream& in, std::ostream& out, std::ostream& err) {
    std::string line;
    getline(in, line);
    while (in && line[0] == '@') {
        out << line << '\n';
        if (&out != &err)
            err << line << '\n';
        getline(in, line);
    }

    std::map<std::string, sam_read*> yet_unpaired;

    while (in) {
        sam_read *r1 = new sam_read(line);
        if (r1->paired()) {
            auto it = yet_unpaired.find(r1->get_name());
            if (it == yet_unpaired.end()) {
                yet_unpaired[r1->get_name()] = r1;
                getline(in, line);
                continue;
            }

            sam_read *r2 = it->second;
            yet_unpaired.erase(it);
            if (r1->first_in_pair() && r2->second_in_pair())
                trim_and_write(*r1, *r2, false, ai, amplicons, out, out, err, err);
            else if (r1->second_in_pair() && r2->first_in_pair())
                trim_and_write(*r2, *r1, false, ai, amplicons, out, out, err, err);
            else {
                trim_and_write(*r1, ai, amplicons, out, err);
                trim_and_write(*r2, ai, amplicons, out, err);
            }
            delete r1;
            delete r2;
        } else {
            trim_and_write(*r1, ai, amplicons, out, err);
            delete r1;
        }
        getline(in, line);
    }

    for (auto const& entry : yet_unpaired) {
        trim_and_write(*entry.second, ai, amplicons, out, err);
    }
}

} // atlas
