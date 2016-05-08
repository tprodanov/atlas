#include "read.h"
#include <sstream>

namespace atlas {

using namespace std;

// READ

read::read()
{}

read::read(string&& name, string&& seq)
    : name_(move(name))
    , seq_(move(seq))
{}

read::~read()
{}

void read::trim(int left, int right, read& other, int o_left, int o_right) {
    trim(left, right);
    other.trim(o_left, o_right);
}

string const& read::get_name() const {
    return name_;
}

string const& read::get_seq() const {
    return seq_;
}

// FASTA

fasta_read::fasta_read(string&& name, string&& seq)
    : read(move(name), move(seq))
{}

void fasta_read::trim(int left, int right) {
    seq_ = seq_.substr(left, right - left);
}

void fasta_read::write_to(ostream& out) const {
    out << name_ << "\n";
    for (size_t i = 1; i <= seq_.size(); ++i) {
        out << seq_[i - 1];
        if (i % 80 == 0) {
            out << "\n";
        }
    }
    out << endl;
}

vector<fasta_read> read_fasta(istream& fasta) {
    string name;
    string seq;
    string current;
    getline(fasta, name);

    vector<fasta_read> reads;

    while (fasta) {
        getline(fasta, current);
        if (current[0] == '>') {
            reads.push_back(fasta_read(std::move(name), std::move(seq)));
            seq = "";
            name = std::move(current);
        } else {
            seq += current;
        }
    }
    reads.push_back(fasta_read(std::move(name), std::move(seq)));
    return reads;
}

// FASTQ

fastq_read::fastq_read()
{}

fastq_read::fastq_read(string&& name, string&& seq, string&& qual)
    : read(move(name), move(seq))
    , qual_(move(qual))
{}

bool fastq_read::from_stream(istream& in, fastq_read& r) {
    getline(in, r.name_);
    if (!in)
        return false;
    getline(in, r.seq_);
    string tmp;
    getline(in, tmp);
    getline(in, r.qual_);
    return true;
}

void fastq_read::trim(int left, int right) {
    seq_ = seq_.substr(left, right - left);
    qual_ = qual_.substr(left, right - left);
}

void fastq_read::write_to(ostream &out) const {
    out << name_ << '\n';
    out << seq_ << '\n';
    out << "+\n";
    out << qual_ << endl;
}

// SAM

sam_read::sam_read()
{}

void from_line(std::string& dest, std::string const& line, size_t& i, size_t& j) {
    j = line.find_first_of('\t', i);
    dest = line.substr(i, j - i);
    i = j + 1;
}

void from_line(int& dest, std::string const& line, size_t& i, size_t& j) {
    j = line.find_first_of('\t', i);
    dest = stoi(line.substr(i, j - i));
    i = j + 1;
}

void sam_read::init(std::string const& line) {
    size_t i = 0;
    size_t j;

    from_line(name_, line, i, j);
    from_line(flag_,  line, i, j);
    from_line(rname_, line, i, j);
    from_line(pos_,   line, i, j);
    from_line(mapq_,  line, i, j);
    from_line(cigar_, line, i, j);
    from_line(rnext_, line, i, j);
    from_line(pnext_, line, i, j);
    from_line(tlen_,  line, i, j);
    from_line(seq_,   line, i, j);
    from_line(qual_,  line, i, j);
    optional_ = line.substr(j);
}

inline bool cigar_length_connected(char symbol) {
    return symbol == 'M' ||
           symbol == 'I' ||
           symbol == 'S' ||
           symbol == '=' ||
           symbol == 'X';
}

std::string trim_cigar(std::string const& cigar, int l_trim_size, int r_trim_size) {
    size_t i = 0;
    int left_number = 0;
    char left_symbol = '\0';
    while (l_trim_size) {
        left_number = 0;

        while ('0' <= cigar[i] && cigar[i]<= '9') {
            left_number = left_number * 10 + (cigar[i] - '0');
            ++i;
        }
        left_symbol = cigar[i];
        ++i;

        if (cigar_length_connected(left_symbol)) {
            int decrease = std::min(l_trim_size, left_number);
            l_trim_size -= decrease;
            left_number -= decrease;
        }
    }

    size_t j = cigar.size() - 1;
    int right_number = 0;
    char right_symbol = '\0';
    while (r_trim_size) {
        right_symbol = cigar[j];
        --j;
        int power = 1;
        right_number = 0;

        while ('0' <= cigar[j] && cigar[j] <= '9') {
            right_number += (cigar[j] - '0') * power;
            power *= 10;
            --j;
        }

        if (cigar_length_connected(right_symbol)) {
            int decrease = std::min(r_trim_size, right_number);
            r_trim_size -= decrease;
            right_number -= decrease;
        }
    }

    std::ostringstream oss;
    if (left_number)
        oss << left_number << left_symbol;
    oss << cigar.substr(i, j + 1 - i);
    if (right_number)
        oss << right_number << right_symbol;
    return oss.str();
}

void sam_read::trim(int left, int right) {
    pos_ += left;

    if (cigar_ != "*")
        cigar_ = trim_cigar(cigar_, left, seq_.size() - right);
    seq_ = seq_.substr(left, right - left);
    qual_ = qual_.substr(left, right - left);
}

void sam_read::trim(int left, int right, read& other, int o_left, int o_right) {
    trim(left, right);
    sam_read& casted = static_cast<sam_read&>(other);
    casted.trim(o_left, o_right);
    pnext_ = casted.pos_;
    casted.pnext_ = pos_;
}

void sam_read::write_to(std::ostream& out) const {
    out << name_ << '\t' << flag_ << '\t' << rname_ << '\t' << pos_ << '\t' << mapq_
        << '\t' << cigar_ << '\t' << rnext_ << '\t' << pnext_ << '\t' << tlen_
        << '\t' << seq_ << '\t' << qual_ << '\t' << optional_ << std::endl;
}

bool sam_read::paired() const {
    return pnext_;
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
    case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M': case 'B': case 'D': case 'H': case 'V':
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
        case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M': case 'B': case 'D': case 'H': case 'V':
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
