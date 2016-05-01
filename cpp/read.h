#pragma once

#include <string>
#include <istream>
#include <ostream>
#include <vector>

namespace atlas {

using u64 = unsigned long long;

struct read {

    read();
    read(std::string const& name, std::string const& seq, std::string const& qual);

    void write_to_fasta(std::ostream& fasta) const;
    void write_to_fastq(std::ostream& fastq) const;

    void trim(int left, int right);

    std::string const& get_name() const;
    std::string const& get_seq() const;
    std::string const& get_qual() const;
    std::size_t size() const noexcept;

    static bool from_fastq(std::istream& fastq, read& r);

private:

    std::string name_;
    std::string seq_;
    std::string qual_;

};

struct unexpected_nt : std::exception {};

std::vector<read> read_fasta(std::istream& fasta);
int get_nt_number(char nt);
std::string reverse_complement(std::string const& seq);
bool seq_kmer(std::string const& seq, size_t k, size_t begin, u64& kmer);

} // namespace atlas
