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

void exit_with_help(std::string const& usage, int exit_code = 64) {
    std::cerr << usage << std::endl;
    exit(exit_code);
}

void exit_with_help(std::string const& error, std::string const& usage, int exit_code = 64) {
    std::cerr << "Error: " << error << "\n\n" << usage << std::endl;
    exit(exit_code);
}

std::string to_upper(std::string const& str) {
    static std::locale loc;
    std::string upper;
    for (char ch : str) {
        upper.push_back(std::toupper(ch, loc));
    }
    return upper;
}

std::istream* get_istream(std::string const& dest, optparse::Values& options,
                          optparse::OptionParser const& parser, std::istream* def = nullptr) {
    if (!options.is_set(dest)) {
        if (def) {
            return def;
        } else {
            exit_with_help(to_upper(dest) + " is required argument", parser.format_help());
        }
    }

    std::istream *stream = new std::ifstream(options[dest]);
    if (!*stream) {
        std::cerr << "Cannot open " << to_upper(dest) << ": " << options[dest] << std::endl;
        exit(64);
    }
    return stream;
}

std::ostream* get_ostream(std::string const& dest, optparse::Values& options,
                          optparse::OptionParser const& parser, std::ostream* def = nullptr) {
    if (!options.is_set(dest)) {
        if (def) {
            return def;
        } else {
            exit_with_help(to_upper(dest) + " is required argument", parser.format_help());
        }
    }

    std::ostream *stream = new std::ofstream(options[dest]);
    if (!*stream) {
        std::cerr << "Cannot write to " << to_upper(dest) << ": " << options[dest] << std::endl;
        exit(64);
    }
    return stream;
}

void build(int argc, char const** argv) {
    optparse::OptionParser parser = optparse::OptionParser()
            .description("Build index");

    parser.add_option("-a", "--amplicons")
            .dest("amplicons")
            .help("amplicons in fasta format [default: stdin]");
    parser.add_option("-i", "--index")
            .dest("index")
            .help("output: amplicon index [default: stdout]");
    parser.add_option("-k")
            .dest("k")
            .set_default("13")
            .help("kmer size [default: 13]");

    std::string name = std::string(argv[0]) + " build";
    argv[1] = name.c_str();
    optparse::Values options = parser.parse_args(argc - 1, argv + 1);

    std::istream *amplicons_stream = get_istream("amplicons", options, parser, &std::cin);
    std::ostream *index_stream = get_ostream("index", options, parser, &std::cout);

    size_t k = (int)options.get("k");
    if (k > 32) {
        std::cerr << "k (" << k << ") should be <= 32" << std::endl;
        exit(64);
    }

    auto amplicons = atlas::create_amplicon_pairs(*amplicons_stream);
    atlas::amplicon_index amplicon_index(amplicons, k);
    amplicon_index.save(*index_stream);

    if (amplicons_stream != &std::cin) {
        delete amplicons_stream;
    }
    if (index_stream != &std::cout) {
        delete index_stream;
    }
}

void trim1(int argc, char const** argv) {
    optparse::OptionParser parser = optparse::OptionParser()
            .description("Trim single reads");

    parser.add_option("-a", "--amplicons")
            .dest("amplicons")
            .help("amplicons in fasta format");
    parser.add_option("-i", "--index")
            .dest("index")
            .help("amplicon index");
    parser.add_option("-f", "--fastq")
            .dest("fastq")
            .help("reads in fastq format [default: stdin]");
    parser.add_option("-o", "--output")
            .dest("output")
            .help("output: trimmed reads [default: stdout]");
    parser.add_option("-e", "--error")
            .dest("error")
            .help("output: unmatched reads [default: stderr]");

    parser.add_option("-t", "--threshold")
            .dest("threshold")
            .set_default(80)
            .help("(optional) alignment score threshold [default: 80]");
    parser.add_option("-m", "--mismatch")
            .dest("mismatch")
            .set_default(-1.4f)
            .help("(optional) mismatch penalty [default: -1.4]");
    parser.add_option("-g", "--gap")
            .dest("gap")
            .set_default(-1.4f)
            .help("(optional) gap penalty [default: -1.4]");

    std::string name = std::string(argv[0]) + " trim1";
    argv[1] = name.c_str();
    optparse::Values options = parser.parse_args(argc - 1, argv + 1);

    std::istream *amplicons_stream = get_istream("amplicons", options, parser);
    std::istream *index_stream = get_istream("index", options, parser);
    std::istream *fastq_stream = get_istream("fastq", options, parser, &std::cin);

    std::ostream *output_stream = get_ostream("output", options, parser, &std::cout);
    std::ostream *error_stream = get_ostream("error", options, parser, &std::cerr);

    atlas::trim_parameters::score_threshold = options.get("threshold");
    atlas::trim_parameters::mismatch_penalty = options.get("mismatch");
    atlas::trim_parameters::gap_penalty = options.get("gap");

    auto amplicons = atlas::create_amplicon_pairs(*amplicons_stream);
    atlas::amplicon_index amplicon_index(*index_stream);

    atlas::trim(amplicon_index, amplicons, *fastq_stream, *output_stream, *error_stream);

    delete amplicons_stream;
    delete index_stream;
    if (fastq_stream != &std::cin)
        delete fastq_stream;
    if (output_stream != &std::cout)
        delete output_stream;
    if (error_stream != &std::cerr)
        delete error_stream;
}


void trim2(int argc, char const** argv) {
    optparse::OptionParser parser = optparse::OptionParser()
            .description("Trim paired reads");

    parser.add_option("-a", "--amplicons")
            .dest("amplicons")
            .help("amplicons in fasta format");
    parser.add_option("-i", "--index")
            .dest("index")
            .help("amplicon index");
    parser.add_option("-1")
            .dest("fastq1");
    parser.add_option("-2")
            .dest("fastq2")
            .help("reads in fastq format");
    parser.add_option("-o", "--output1")
            .dest("output1");
    parser.add_option("-O", "--output2")
            .dest("output2")
            .help("output: trimmed reads");
    parser.add_option("-e", "--error1")
            .dest("error1");
    parser.add_option("-E", "--error2")
            .dest("error2")
            .help("output: unmatched reads");

    parser.add_option("-t", "--threshold")
            .dest("threshold")
            .set_default(80)
            .help("(optional) alignment score threshold [default: 80]");
    parser.add_option("-m", "--mismatch")
            .dest("mismatch")
            .set_default(-1.4f)
            .help("(optional) mismatch penalty [default: -1.4]");
    parser.add_option("-g", "--gap")
            .dest("gap")
            .set_default(-1.4f)
            .help("(optional) gap penalty [default: -1.4]");

    std::string name = std::string(argv[0]) + " trim2";
    argv[1] = name.c_str();
    optparse::Values options = parser.parse_args(argc - 1, argv + 1);

    std::istream *amplicons_stream = get_istream("amplicons", options, parser);
    std::istream *index_stream = get_istream("index", options, parser);
    std::istream *fastq1_stream = get_istream("fastq1", options, parser);
    std::istream *fastq2_stream = get_istream("fastq2", options, parser);

    std::ostream *output1_stream = get_ostream("output1", options, parser);
    std::ostream *output2_stream = get_ostream("output2", options, parser);
    std::ostream *error1_stream = get_ostream("error1", options, parser);
    std::ostream *error2_stream = get_ostream("error2", options, parser);

    atlas::trim_parameters::score_threshold = options.get("threshold");
    atlas::trim_parameters::mismatch_penalty = options.get("mismatch");
    atlas::trim_parameters::gap_penalty = options.get("gap");

    auto amplicons = atlas::create_amplicon_pairs(*amplicons_stream);
    atlas::amplicon_index amplicon_index(*index_stream);

    atlas::trim(amplicon_index, amplicons,
                *fastq1_stream, *fastq2_stream,
                *output1_stream, *output2_stream,
                *error1_stream, *error2_stream);

    delete amplicons_stream;
    delete index_stream;
    delete fastq1_stream;
    delete fastq2_stream;
    delete output1_stream;
    delete output2_stream;
    delete error1_stream;
    delete error2_stream;
}

int main(int argc, char const** argv) {
    std::string name = optparse::basename(argv[0]);
    std::string usage =
    R"(Amplicon Trimmer.

Usage: )" + name + R"( amplicon_trimmer <command> [<args>...]

Possible commands are:
    build   build aplicon index
    trim1   trim single reads
    help    show this message
    )";

    if (argc < 2) {
        exit_with_help(usage, 0);
    }

    std::string command = argv[1];

    if (command == "build") {
        build(argc, argv);
    } else if (command == "trim1") {
        trim1(argc, argv);
    } else if (command == "trim2") {
        trim2(argc, argv);
    } else if (command == "help") {
        exit_with_help(usage, 0);
    } else {
        exit_with_help(usage);
    }

    return 0;
}
