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

void read::trim(size_t left, size_t right, read& other, size_t o_left, size_t o_right) {
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

void fasta_read::trim(size_t left, size_t right) {
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

void fastq_read::trim(size_t left, size_t right) {
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

sam_read::sam_read(std::string const& line) {
    init(line);
}

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
    trimmed_ = false;
    size_t i = 0;
    size_t j;

    from_line(name_, line, i, j);
    from_line(flag_,  line, i, j);
    from_line(rname_, line, i, j);
    from_line(pos_,   line, i, j);
    from_line(mapq_,  line, i, j);
    std::string cigar;
    from_line(cigar, line, i, j);
    construct_cigar(cigar);
    from_line(rnext_, line, i, j);
    from_line(pnext_, line, i, j);
    from_line(tlen_,  line, i, j);
    from_line(seq_,   line, i, j);
    from_line(qual_,  line, i, j);
    optional_ = line.substr(j + 1);
}

void sam_read::construct_cigar(std::string const& cigar) {
    cigar_.clear();
    if (cigar == "*")
        return;

    size_t current_number = 0;
    for (char ch : cigar) {
        if (ch >= '0' && ch <= '9')
            current_number = current_number * 10 + (ch - '0');
        else {
            cigar_.push_back(std::make_pair(current_number, ch));
            current_number = 0;
        }
    }
    clip_cigar();
}

inline bool cigar_length_connected(char symbol) {
    return symbol == 'M' ||
           symbol == 'I' ||
           symbol == 'S' ||
           symbol == '=' ||
           symbol == 'X';
}


void sam_read::clip_cigar() {
    while (!cigar_.empty() && !cigar_length_connected(cigar_.front().second))
        cigar_.pop_front();
    while (!cigar_.empty() && !cigar_length_connected(cigar_.back().second))
        cigar_.pop_back();
}


void sam_read::trim_cigar(size_t l_trim_size, size_t r_trim_size,
                          size_t& l_mapped_trimmed, size_t& r_mapped_trimmed) {
    l_mapped_trimmed = 0;
    while (l_trim_size) {
        if (cigar_length_connected(cigar_.front().second)) {
            size_t decrease = std::min(l_trim_size, cigar_.front().first);
            l_trim_size -= decrease;
            cigar_.front().first -= decrease;

            if (cigar_.front().second != 'S')
                l_mapped_trimmed += decrease;
            if (cigar_.front().first == 0)
                cigar_.pop_front();
        } else
            cigar_.pop_front();
    }

    r_mapped_trimmed = 0;
    while (r_trim_size) {
        if (cigar_length_connected(cigar_.back().second)) {
            size_t decrease = std::min(r_trim_size, cigar_.back().first);
            r_trim_size -= decrease;
            cigar_.back().first -= decrease;

            if (cigar_.back().second != 'S')
                r_mapped_trimmed += decrease;
            if (cigar_.back().first == 0)
                cigar_.pop_back();
        } else
            cigar_.pop_back();
    }
    clip_cigar();
}

void sam_read::trim(size_t left, size_t right, size_t& l_mapped_trimmed, size_t& r_mapped_trimmed) {
    trimmed_ = true;
    if (!cigar_.empty())
        trim_cigar(left, seq_.size() - right, l_mapped_trimmed, r_mapped_trimmed);
    pos_ += l_mapped_trimmed;
    seq_ = seq_.substr(left, right - left);
    qual_ = qual_.substr(left, right - left);
}

void sam_read::trim(size_t left, size_t right) {
    size_t l_mapped_trimmed;
    size_t r_mapped_trimmed;
    trim(left, right, l_mapped_trimmed, r_mapped_trimmed);
}

void sam_read::update_tlen(sam_read& s1, int l1, int r1, sam_read& s2, int l2, int r2) {
    int change1 = s1.is_reverse_complement() ? r1 : -l1;
    int change2 = s2.is_reverse_complement() ? -r2 : l2;
    s1.tlen_ += change1 + change2;
    s2.tlen_ -= change1 + change2;
}


void sam_read::trim(size_t left, size_t right, read& other, size_t o_left, size_t o_right) {
    size_t l1_mapped_trimmed;
    size_t r1_mapped_trimmed;
    trim(left, right, l1_mapped_trimmed, r1_mapped_trimmed);

    sam_read& casted = static_cast<sam_read&>(other);
    size_t l2_mapped_trimmed;
    size_t r2_mapped_trimmed;
    casted.trim(o_left, o_right, l2_mapped_trimmed, r2_mapped_trimmed);

    if (tlen_ > 0)
        update_tlen(*this, l1_mapped_trimmed, r1_mapped_trimmed, casted, l2_mapped_trimmed, r2_mapped_trimmed);
    else if (tlen_ < 0)
        update_tlen(casted, l2_mapped_trimmed, r2_mapped_trimmed, *this, l1_mapped_trimmed, r1_mapped_trimmed);

    pnext_ = casted.pos_;
    casted.pnext_ = pos_;
}

bool delete_flag(char letter1, char letter2) {
    if (letter1 == 'C')
        return letter2 == 'C' || letter2 == 'P';
    if (letter1 == 'H')
        return letter2 >= '0' && letter2 <= '2';
    if (letter1 == 'M')
        return letter2 == 'C' || letter2 == 'D';
    if (letter1 == 'N')
        return letter2 == 'M';
    if (letter1 == 'O')
        return letter2 == 'P' || letter2 == 'Q';

    if (letter2 == '2')
        return letter1 == 'E' || letter1 == 'Q' || letter1 == 'R';
    return false;
}

void sam_read::write_optional(std::ostream& out) const {
    size_t i = 0;
    size_t j = 0;
    while (j < optional_.size()) {
        std::string flag;
        from_line(flag, optional_, i, j);
        if (!delete_flag(flag[0], flag[1]))
            out << '\t' << flag;
    }
}

void sam_read::write_to(std::ostream& out) const {
    out << name_ << '\t' << flag_ << '\t' << rname_ << '\t' << pos_ << '\t' << mapq_
        << '\t';
    for (auto const& entry : cigar_)
        out << entry.first << entry.second;
    if (cigar_.empty())
        out << '*';
    out << '\t' << rnext_ << '\t' << pnext_ << '\t' << tlen_
        << '\t' << seq_ << '\t' << qual_;
    if (trimmed_)
        write_optional(out);
    else
        out << optional_;
    out << std::endl;
}

bool sam_read::paired() const {
    return (flag_ & 0x1) &&         // Has pair
            (flag_ & 0x800) == 0 && // Not supplementary
            pnext_;
}

bool sam_read::first_in_pair() const {
    return flag_ & 0x40;
}

bool sam_read::second_in_pair() const {
    return flag_ & 0x80;
}

bool sam_read::is_reverse_complement() const {
    return flag_ & 0x10;
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

void recursive_n_filler(std::vector<u64> const& n_multipliers, std::vector<u64>& kmers, u64 kmer, size_t pos = 0) {
    if (pos >= n_multipliers.size()) {
        kmers.push_back(kmer);
        return;
    }
    for (u64 i = 0; i < 4; ++i)
        recursive_n_filler(n_multipliers, kmers, kmer + n_multipliers[pos] * i, pos + 1);
}

std::vector<u64> seq_kmer(std::string const& seq, size_t k, size_t begin, size_t d_max) {
    std::vector<u64> kmers;
    std::vector<u64> n_multipliers;

    u64 base_kmer = 0;
    for (size_t i = 0; i < k; ++i) {
        if (seq[i + begin] == 'N') {
            if (n_multipliers.size() >= d_max)
                return kmers;
            n_multipliers.push_back(1 << 2 * i);
            base_kmer *= 4;
        } else
            base_kmer = 4 * base_kmer + get_nt_number(seq[i + begin]);
    }

    if (n_multipliers.empty())
        kmers.push_back(base_kmer);
    else
        recursive_n_filler(n_multipliers, kmers, base_kmer);
    return kmers;
}

} // namespace atlas
