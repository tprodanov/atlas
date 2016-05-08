#pragma once

#include <string>
#include <set>

#include "amplicon_index.h"

namespace atlas {

struct trim_parameters {

    static float score_threshold_ratio;
    static float gap_penalty;
    static float mismatch_penalty;

};

enum direction {
    stop = 0,
    up,
    diagonal,
    left
};

size_t find_start(direction* prev, size_t coord, size_t n);
float trim_to_amplicon(std::string const& seq, std::string const& amplicon, size_t& left, size_t& right);

std::set<int> possible_amplicons(std::string const& seq,
                                 amplicon_index const& ai);
bool trim_read(read& r, amplicon_index const& ai, amplicon_pairs const& amplicons);
bool trim_reads(read& r1, read& r2, amplicon_index const& ai, amplicon_pairs const& amplicons);

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq, std::ostream& output, std::ostream& error);
void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq1,  std::istream& fastq2,
          std::ostream& output1, std::ostream& output2,
          std::ostream& error1,  std::ostream& error2);

void trim_sam(amplicon_index const& ai, amplicon_pairs const& amplicons,
              std::istream& in, std::ostream& out, std::ostream& err);

} // namespace atlas
