#include "argparse/option_parser.h"

#include <iostream>
#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <fstream>
#include <locale>

#include "amplicon_index.h"
#include "trim.h"
#include "parse_args.h"

void trim_fq1(atlas::parse_args const& parser, atlas::amplicon_pairs const& amplicons,
              atlas::amplicon_index const& amplicon_index) {
    std::istream *fastq_stream = parser.get_istream("f", "<in.fq>", 1, &std::cin);
    std::ostream *output_stream = parser.get_ostream("o", "<out.fq>", 1, &std::cout);
    std::ostream *error_stream = parser.get_ostream("e", "<err.fq>", 1, output_stream);

    atlas::trim(amplicon_index, amplicons, *fastq_stream, *output_stream, *error_stream);

    if (fastq_stream != &std::cin)
        delete fastq_stream;
    if (error_stream != output_stream)
        delete error_stream;
    if (output_stream != &std::cout)
        delete output_stream;
}

void trim_fq2(atlas::parse_args const& parser, atlas::amplicon_pairs const& amplicons,
              atlas::amplicon_index const& amplicon_index) {
    std::istream *fastq1_stream = parser.get_istream("f", "<in1.fq>", 1);
    std::istream *fastq2_stream = parser.get_istream("f", "<in2.fq>", 2);
    std::ostream *output1_stream = parser.get_ostream("o", "<out1.fq>", 1);
    std::ostream *output2_stream = parser.get_ostream("o", "<out2.fq>", 2);
    std::ostream *error1_stream = parser.get_ostream("e", "<err1.fq>", 1, output1_stream);
    std::ostream *error2_stream = parser.get_ostream("e", "<err2.fq>", 2, output2_stream);

    atlas::trim_fastq(amplicon_index, amplicons,
                *fastq1_stream, *fastq2_stream,
                *output1_stream, *output2_stream,
                *error1_stream, *error2_stream);

    delete fastq1_stream;
    delete fastq2_stream;
    if (error1_stream != output1_stream)
        delete error1_stream;
    if (error2_stream != output2_stream)
        delete error2_stream;
    delete output1_stream;
    delete output2_stream;
}

void trim_sam(atlas::parse_args const& parser, atlas::amplicon_pairs const& amplicons,
              atlas::amplicon_index const& amplicon_index) {
    std::istream *in_stream = parser.get_istream("s", "<in.sam>", 1, &std::cin);
    std::ostream *output_stream = parser.get_ostream("o", "<out.sam>", 1, &std::cout);
    std::ostream *error_stream = parser.get_ostream("e", "<err.sam>", 1, output_stream);

    atlas::trim_sam(amplicon_index, amplicons, *in_stream, *output_stream, *error_stream);

    if (in_stream != &std::cin)
        delete in_stream;
    if (error_stream != output_stream)
        delete error_stream;
    if (output_stream != &std::cout)
        delete output_stream;
}


void trim(int argc, char const** argv) {
    std::string name = argv[0];
    std::string usage =
    R"(Amplicon Trimmer.

Usage:
    )" + name + R"( trim [options] -a <ampl.fa> -f [<in.fq>] [-o <out.fq>] [-e <err.fq>]
    )" + name + R"( trim [options] -a <ampl.fa> -f <in1.fq> <in2.fq> -o <out1.fq> <out2.fq> [-e <err1.fq> <err2.fq>]
    )" + name + R"( trim [options] -a <ampl.fa> -s [<in.sam>] [-o <out.sam>] [-e <err.sam>]

    .fa indicates fasta format
    .fq indicates fastq format
    .sam indicates sam format

    <ampl.fa>   Amplicons
    <in.*>      File with unpaired reads (default: stdin)
    <in1.*>     File with #1 mates
    <in2.*>     File with #2 mates

    <out.*>     Trimmed unpaired reads (default: stdout)
    <out1.*>    Trimmed #1 mates
    <out2.*>    Trimmed #2 mates
    <err.*>     Unmatched unpaired reads (default: <out.*>)
    <err1.*>    Unmatched #1 mates (default: <out1.*>)
    <err2.*>    Unmatched #2 mates (default: <out2.*>)

Options:
    -k INT      K-mer size (default: 13)
    -t FLOAT    Alignment score ratio threshold (default: 0.5)
    -m FLOAT    Mismatch penalty (default: -1.4)
    -g FLOAT    Gap penalty (default: -1.4)
    -d INT      Max number of degenerate nucleotides N in a k-mer (default: 2)
    -b          Trim paired reads from both sides

Other:
    -h/--help       Show this message
)";

    atlas::parse_args parser(usage);
    parser.specify_min_nargs("m", 1);
    parser.specify_min_nargs("g", 1);
    parser.specify_min_nargs("t", 1);

    parser.parse(argc, argv);

    size_t k = parser.get_value<size_t>("k", "k", 1, 13);
    if (k > 32) {
        std::cerr << "k (" << k << ") should be <= 32" << std::endl;
        exit(64);
    }
    atlas::trim_parameters::score_threshold_ratio = parser.get_value<float>("t", "-t THRESHOLD", 1, 0.5);
    if (atlas::trim_parameters::score_threshold_ratio < 0 ||
            1 < atlas::trim_parameters::score_threshold_ratio) {
        std::cerr << "threshold ratio (" << atlas::trim_parameters::score_threshold_ratio
                  << ") should be in [0, 1]" << std::endl;
        exit(64);
    }
    atlas::trim_parameters::mismatch_penalty = parser.get_value<float>("m", "-m MISMATCH", 1, -1.4);
    atlas::trim_parameters::gap_penalty = parser.get_value<float>("g", "-g GAP", 1, -1.4);
    atlas::trim_parameters::degenerate_nucleotides_count = parser.get_value<size_t>("d", "-d INT", 1, 2);
    atlas::trim_parameters::paired_local = parser.flag_up("b");

    std::istream *amplicons_stream = parser.get_istream("a", "<ampl.fa>", 1);
    auto amplicons = atlas::create_amplicon_pairs(*amplicons_stream);
    delete amplicons_stream;
    atlas::amplicon_index amplicon_index(amplicons, k);

    if (parser.flag_up("f")) {
        if (parser.nargs("f") <= 1) {
            trim_fq1(parser, amplicons, amplicon_index);
        } else {
            trim_fq2(parser, amplicons, amplicon_index);
        }
    } else if (parser.flag_up("s")) {
        trim_sam(parser, amplicons, amplicon_index);
    } else {
        parser.exit_with_help("Either -f or -s flag should be specified");
    }
}

int main(int argc, char const** argv) {
    std::string name = atlas::base_name(argv[0]);
    argv[0] = name.c_str();
    std::string usage =
    R"(Amplicon Trimmer.

Usage: )" + name + R"( amplicon_trimmer <command> [<args>...]

Possible commands are:
    trim    trim reads
    help    show this message
    )";

    if (argc < 2) {
        std::cout << usage << std::endl;
        exit(0);
    }

    std::string command = argv[1];

    if (command == "trim") {
        trim(argc, argv);
    } else {
        std::cout << usage << std::endl;
        exit(0);
    }

    return 0;
}
