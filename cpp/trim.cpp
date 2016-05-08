#include "trim.h"

namespace atlas {

float trim_parameters::score_threshold_ratio = 0.5;
float trim_parameters::gap_penalty = -1.4;
float trim_parameters::mismatch_penalty = -1.4;

size_t find_start(direction* prev, size_t coord, size_t n) {
    while (true) {
        switch (prev[coord]) {
        case stop:
            return coord;
        case left:
            --coord;
            break;
        case up:
            coord -= n + 1;
            break;
        case diagonal:
            coord -= n + 2;
            break;
        }
    }
}

float local_align(std::string const& seq, std::string const& amplicon, size_t& l, size_t& r) {
    size_t n = seq.size();
    size_t m = amplicon.size();
    float *h = new float[(n + 1) * (m + 1)];
    direction *prev = new direction[(n + 1) * (m + 1)];

    float score = -1e6;
    size_t best_coord = 0;

    for (size_t j = 0; j <= n; ++j) {
        h[j] = 0;
        prev[j] = direction::stop;
    }

    for (size_t i = 1; i <= m; ++i) {
        h[i * (n + 1)] = 0;
        prev[i * (n + 1)] = direction::stop;

        for (size_t j = 1; j <= n; ++j) {
            h[i * (n + 1) + j] = h[i * (n + 1) + j - 1] + trim_parameters::gap_penalty;
            prev[i * (n + 1) + j] = direction::left;

            if (h[(i - 1) * (n + 1) + j] + trim_parameters::gap_penalty > h[i * (n + 1) + j]) {
                h[i * (n + 1) + j] = h[(i - 1) * (n + 1) + j] + trim_parameters::gap_penalty;
                prev[i * (n + 1) + j] = direction::up;
            }

            if (h[(i - 1) * (n + 1) + j - 1]
                    + (seq[j - 1] == amplicon[i - 1] ? 1.0f : trim_parameters::mismatch_penalty)
                    > h[i * (n + 1) + j]) {
                h[i * (n + 1) + j] = h[(i - 1) * (n + 1) + j - 1]
                        + (seq[j - 1] == amplicon[i - 1] ? 1.0f : trim_parameters::mismatch_penalty);
                prev[i * (n + 1) + j] = direction::diagonal;
            }

            if (h[i * (n + 1) + j] < 0) {
                h[i * (n + 1) + j] = 0;
                prev[i * (n + 1) + j] = direction::stop;
            }

            if (h[i * (n + 1) + j] > score) {
                score = h[i * (n + 1) + j];
                best_coord = i * (n + 1) + j;
            }
        }
    }

    l = find_start(prev, best_coord, n) % (n + 1);
    r = best_coord % (n + 1);

    delete[] prev;
    delete[] h;
    return score;
}

std::set<int> possible_amplicons(std::string const& seq, amplicon_index const& ai) {
    size_t n = seq.size();
    size_t k = ai.get_k();
    if (n < k) {
        return std::set<int>();
    }

    std::vector<size_t> kmer_pos = { 0, n - k };
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
        if (possible_amplicons1.find(-ix) != possible_amplicons1.end()) {
            size_t left1;
            size_t right1;
            float score1 = ix > 0
                ? local_align(seq1, amplicons[ix - 1].second,   left1, right1)
                : local_align(seq1, amplicons[-ix - 1].first, left1, right1);

            size_t left2;
            size_t right2;
            float score2 = ix > 0
                ? local_align(seq2, amplicons[ix - 1].first, left2, right2)
                : local_align(seq2, amplicons[-ix - 1].second, left2, right2);

            if (score1 + score2 > best_score) {
                best_score = score1 + score2;
                begin1 = left1;
                end1 = right1;

                begin2 = left2;
                end2 = right2;
            }
        }
        if (possible_amplicons1.find(ix) != possible_amplicons1.end()) {
            size_t left1;
            size_t right1;
            float score1 = ix > 0
                ? local_align(seq1, amplicons[ix - 1].first, left1, right1)
                : local_align(seq1, amplicons[-ix - 1].second, left1, right1);

            size_t left2;
            size_t right2;
            float score2 = ix > 0
                ? local_align(seq2, amplicons[ix - 1].first, left2, right2)
                : local_align(seq2, amplicons[-ix - 1].second, left2, right2);

            if (score1 + score2 > best_score) {
                best_score = score1 + score2;
                begin1 = left1;
                end1 = right1;

                begin2 = left2;
                end2 = right2;
            }
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

void trim_and_write(read& r1, read& r2, amplicon_index const& ai, amplicon_pairs const& amplicons,
                    std::ostream& out1, std::ostream& out2, std::ostream& err1, std::ostream& err2) {
    try {
        if (trim_reads(r1, r2, ai, amplicons)) {
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

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq1,  std::istream& fastq2,
          std::ostream& output1, std::ostream& output2,
          std::ostream& error1,  std::ostream& error2) {
    fastq_read r1;
    fastq_read r2;

    while (fastq_read::from_stream(fastq1, r1) && fastq_read::from_stream(fastq2, r2)) {
        trim_and_write(r1, r2, ai, amplicons, output1, output2, error1, error2);
    }
}

void trim_sam(amplicon_index const& ai, amplicon_pairs const& amplicons,
              std::istream& in, std::ostream& out, std::ostream& err) {
    sam_read r1;
    sam_read r2;

    std::string line;
    getline(in, line);
    while (in && line[0] == '@') {
        out << line << '\n';
        if (&out != &err)
            err << line << '\n';
        getline(in, line);
    }

    while (in) {
        r1.init(line);
        if (r1.paired()) {
            getline(in, line);
            r2.init(line);
            trim_and_write(r1, r2, ai, amplicons, out, out, err, err);
        } else {
            trim_and_write(r1, ai, amplicons, out, err);
        }
        getline(in, line);
    }
}

} // atlas
