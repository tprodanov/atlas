#include "trim.h"

namespace atlas {

int trim_parameters::score_threshold = 80;
float trim_parameters::gap_penalty = -1.4;
float trim_parameters::mismatch_penalty = -1.4;

size_t align_from_the_left(std::string const& seq, std::string const& amplicon, float& score) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    score = -1e6;
    size_t best_read_end = 0;

    for (size_t i = 0; i <= m; ++i) {
        h[i * (n + 1)] = 0;
        for (size_t j = 1; j <= n; ++j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j - 1] + trim_parameters::gap_penalty;

            if (i) {
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j], h[(i - 1) * (n + 1) + j]
                        + trim_parameters::gap_penalty);
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j], h[(i - 1) * (n + 1) + j - 1]
                        + (seq[i - 1] == amplicon[j - 1] ? 1.0f : trim_parameters::mismatch_penalty));
            }

            if (h[i * (n + 1) + j] > score) {
                score = h[i * (n + 1) + j];
                best_read_end = j;
            }
        }
    }
    delete[] h;
    return best_read_end;
}

size_t align_from_the_right(std::string const& seq, std::string const& amplicon, float& score) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    score = -1e6;
    size_t best_read_begin = n;

    for (size_t i = m; i <= m; --i) {
        h[i * (n + 1) + n] = 0;
        for (size_t j = n - 1; j < n; --j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j + 1] + trim_parameters::gap_penalty;

            if (i < m) {
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j], h[(i + 1) * (n + 1) + j]
                        + trim_parameters::gap_penalty);
                h[i * (n + 1) + j] = std::max(h[i * (n + 1) + j], h[(i + 1) * (n + 1) + j + 1]
                        + (seq[i] == amplicon[j] ? 1.0f : trim_parameters::mismatch_penalty));
            }

            if (h[i * (n + 1) + j] > score) {
                score = h[i * (n + 1) + j];
                best_read_begin = j;
            }
        }
    }
    delete[] h;
    return best_read_begin;
}

float trim_to_amplicon(std::string const& seq, std::string const& amplicon, size_t& left, size_t& right) {
    size_t n = seq.size();

    float l_score;
    size_t read_end = align_from_the_left(seq, amplicon, l_score);

    float r_score;
    size_t read_begin = align_from_the_right(seq, amplicon, r_score);

    if (l_score < r_score) {
        left = read_begin;
        right = n;
        return r_score;
    } else {
        left = 0;
        right = read_end;
        return l_score;
    }
}

std::set<int> possible_amplicons(std::string const& seq, amplicon_index const& ai) {
    size_t n = seq.size();
    size_t k = ai.get_k();
    if (n <= k) {
        return std::set<int>();
    }

    std::vector<size_t> kmer_pos = { 0, k, n - k };
    if (n >= 2 * k) {
        kmer_pos.push_back(k);
        kmer_pos.push_back(n / 2 - k / 2);
        kmer_pos.push_back(n - 2 * k);
    }

    std::set<int> result;
    for (size_t pos : kmer_pos) {
        u64 kmer;
        if (seq_kmer(seq, k, pos, kmer)) {
            if (!ai.contains(kmer)) {
                continue;
            }
            for (int amplicon : ai.get_all_amplicons(kmer)) {
                result.insert(amplicon);
            }
        }
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
            ? trim_to_amplicon(seq, amplicons[ix - 1].first, left, right)
            : trim_to_amplicon(seq, amplicons[-ix - 1].second, left, right);
        if (score > best_score) {
            best_score = score;
            best_begin = left;
            best_end = right;
        }
    }

    if (best_score >= trim_parameters::score_threshold) {
        r.trim(best_begin, best_end);
        return true;
    }
    return false;
}

bool trim_reads(read& r1, read& r2, amplicon_index const& ai, amplicon_pairs const& amplicons) {
    std::string const& seq1 = r1.get_seq();
    std::string const& seq2 = r2.get_seq();
    float best_score = -1;

    size_t begin1 = 0;
    size_t end1 = seq1.size();

    size_t begin2 = 0;
    size_t end2 = seq2.size();

    std::set<int> possible_amplicons1 = possible_amplicons(seq1, ai);

    for (int ix : possible_amplicons(seq2, ai)) {
        if (possible_amplicons1.find(-ix) == possible_amplicons1.end()) {
            continue;
        }

        size_t left1;
        size_t right1;
        float score1 = ix > 0
            ? trim_to_amplicon(seq1, amplicons[ix - 1].second,   left1, right1)
            : trim_to_amplicon(seq1, amplicons[-ix - 1].first, left1, right1);

        size_t left2;
        size_t right2;
        float score2 = ix > 0
            ? trim_to_amplicon(seq2, amplicons[ix - 1].first, left2, right2)
            : trim_to_amplicon(seq2, amplicons[-ix - 1].second, left2, right2);

        if (score1 + score2 > best_score) {
            best_score = score1 + score2;
            begin1 = left1;
            end1 = right1;

            begin2 = left2;
            end2 = right2;
        }
    }

    if (best_score >= 2 * trim_parameters::score_threshold) {
        r1.trim(begin1, end1);
        r2.trim(begin2, end2);
        return true;
    }
    return false;

}

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq, std::ostream& output, std::ostream& error) {
    read r;
    while (read::from_fastq(fastq, r)) {
        if (trim_read(r, ai, amplicons)) {
            r.write_to_fastq(output);
        } else {
            r.write_to_fastq(error);
        }
    }
}

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq1,  std::istream& fastq2,
          std::ostream& output1, std::ostream& output2,
          std::ostream& error1,  std::ostream& error2) {
    read r1;
    read r2;

    while (read::from_fastq(fastq1, r1) && read::from_fastq(fastq2, r2)) {
        if (trim_reads(r1, r2, ai, amplicons)) {
            r1.write_to_fastq(output1);
            r2.write_to_fastq(output2);
        } else {
            r1.write_to_fastq(error1);
            r2.write_to_fastq(error2);
        }
    }
}

} // atlas
