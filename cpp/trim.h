#pragma once

#include <string>
#include <set>

#include "amplicon_index.h"

namespace atlas {

struct trim_parameters {

    static float score_threshold_ratio;
    static float gap_penalty;
    static float mismatch_penalty;
    static size_t degenerate_nucleotides_count;
    static bool paired_local;

};

void trim(amplicon_index const& ai, amplicon_pairs const& amplicons,
          std::istream& fastq, std::ostream& output, std::ostream& error);
void trim_fastq(amplicon_index const& ai, amplicon_pairs const& amplicons,
                std::istream& fastq1,  std::istream& fastq2,
                std::ostream& output1, std::ostream& output2,
                std::ostream& error1,  std::ostream& error2);

void trim_sam(amplicon_index const& ai, amplicon_pairs const& amplicons,
              std::istream& in, std::ostream& out, std::ostream& err);

} // namespace atlas
