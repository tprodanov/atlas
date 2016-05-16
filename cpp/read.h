#pragma once

#include <string>
#include <istream>
#include <ostream>
#include <vector>
#include <deque>

namespace atlas {

using u64 = unsigned long long;

struct read {

    read();
    read(std::string&& name, std::string&& seq);

    virtual ~read();
    virtual void trim(size_t left, size_t right) = 0;
    virtual void trim(size_t left, size_t right, read& other, size_t o_left, size_t o_right);

    virtual void write_to(std::ostream& out) const = 0;
    virtual std::string const& get_name() const;
    virtual std::string const& get_seq() const;

protected:

    std::string name_;
    std::string seq_;

};

struct fasta_read
    : public read {

    fasta_read(std::string&& name, std::string&& seq);

    virtual void trim(size_t left, size_t right) override;
    virtual void write_to(std::ostream& out) const override;

};

std::vector<fasta_read> read_fasta(std::istream& fasta);

struct fastq_read
    : public read {

    fastq_read();
    fastq_read(std::string&& name, std::string&& seq, std::string&& qual);
    static bool from_stream(std::istream& in, fastq_read& r);

    virtual void trim(size_t left, size_t right) override;
    virtual void write_to(std::ostream& out) const override;

    std::string const& get_qual() const;

private:

    std::string qual_;

};

struct sam_read
    : public read {

    sam_read();
    sam_read(std::string const& line);
    void init(std::string const& line);

    virtual void trim(size_t left, size_t right) override;
    virtual void trim(size_t left, size_t right, read& other, size_t o_left, size_t o_right) override;
    virtual void write_to(std::ostream& out) const override;

    bool paired() const;
    bool first_in_pair() const;
    bool second_in_pair() const;

private:

    void construct_cigar(std::string const& cigar);
    void clip_cigar();
    void trim_cigar(size_t l_trim_size, size_t r_trim_size, size_t& l_mapped_trimmed, size_t& r_mapped_trimmed);
    bool is_reverse_complement() const;
    void trim(size_t left, size_t right, size_t& l_mapped_trimmed, size_t& r_mapped_trimmed);
    static void update_tlen(sam_read& s1, int l1, int r1, sam_read& s2, int l2, int r2);
    void write_optional(std::ostream& out) const;

    int         flag_;
    std::string rname_;
    int         pos_;
    int         mapq_;
    std::deque<std::pair<size_t, char>> cigar_;
    std::string rnext_;
    int         pnext_;
    int         tlen_;
    std::string qual_;
    std::string optional_;
    bool trimmed_;

};

struct unexpected_nt : std::exception {};


int get_nt_number(char nt);
std::string reverse_complement(std::string const& seq);
std::vector<u64> seq_kmer(std::string const& seq, size_t k, size_t begin, size_t d_max);

} // namespace atlas
