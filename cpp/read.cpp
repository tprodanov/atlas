#include "read.h"

namespace atlas {

using namespace std;

read::read()
{}

read::read(string const& name, string const& seq, string const& qual)
    : name_(name)
    , seq_(seq)
    , qual_(qual)
{}

bool read::from_fastq(istream& fastq, read& r) {
    getline(fastq, r.name_);
    if (!fastq)
        return false;
    r.name_ = r.name_.substr(1);
    getline(fastq, r.seq_);
    string tmp;
    getline(fastq, tmp);
    getline(fastq, r.qual_);
    return true;
}

vector<read> read_fasta(istream& fasta) {
    string name;
    string seq;
    string current;
    getline(fasta, name);

    vector<read> reads;

    while (fasta) {
        getline(fasta, current);
        if (current[0] == '>') {
            reads.push_back(read(name.substr(1), std::move(seq), ""));
            seq = "";
            name = std::move(current);
        } else {
            seq += current;
        }
    }
    reads.push_back(read(name.substr(1), std::move(seq), ""));
    return reads;
}

void read::write_to_fasta(ostream& fasta) const {
    fasta << ">" << name_ << "\n";
    for (size_t i = 1; i <= seq_.size(); ++i) {
        fasta << seq_[i - 1];
        if (i % 80 == 0) {
            fasta << "\n";
        }
    }
    fasta << endl;
}

void read::write_to_fastq(std::ostream& fastq) const {
    fastq << "@" << name_ << "\n";
    fastq << seq_ << "\n";
    fastq << "+\n";
    fastq << qual_ << endl;
}

void read::trim(int left, int right) {
    seq_ = seq_.substr(left, right);
    if (qual_.size()) {
        qual_ = qual_.substr(left, right);
    }
}

string const& read::get_name() const {
    return name_;
}

string const& read::get_seq() const {
    return seq_;
}

string const& read::get_qual() const {
    return qual_;
}

size_t read::size() const noexcept {
    return seq_.size();
}

int get_nt_number(char nt) {
    switch (nt) {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    case 'N':
        return -1;
    }
    throw unexpected_nt();
}

std::string reverse_complement(std::string const& seq) {
    std::string rc;
    for (size_t i = seq.size(); i > 0; --i) {
        switch (seq[i - 1]) {
        case 'A':
            rc.push_back('T');
            break;
        case 'C':
            rc.push_back('G');
            break;
        case 'G':
            rc.push_back('C');
            break;
        case 'T':
            rc.push_back('A');
            break;
        case 'N':
            rc.push_back('N');
            break;
        default:
            throw unexpected_nt();
        }
    }
    return rc;
}

bool seq_kmer(std::string const& seq, size_t k, size_t begin, u64& kmer) {
    kmer = 0;
    for (size_t i = 0; i < k; ++i) {
        if (seq[i + begin] == 'N')
            return false;
        kmer = 4 * kmer + get_nt_number(seq[i + begin]);
    }
    return kmer;
}

} // namespace atlas
