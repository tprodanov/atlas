#pragma once

#include <string>
#include <istream>
#include <ostream>
#include <vector>

namespace atlas {

using u64 = unsigned long long;

struct read {

    read();
    read(std::string&& name, std::string&& seq);

    virtual ~read();
    virtual void trim(int left, int right) = 0;
    virtual void trim(int left, int right, read& other, int o_left, int o_right);

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

    virtual void trim(int left, int right) override;
    virtual void write_to(std::ostream& out) const override;

};

std::vector<fasta_read> read_fasta(std::istream& fasta);

struct fastq_read
    : public read {

    fastq_read();
    fastq_read(std::string&& name, std::string&& seq, std::string&& qual);
    static bool from_stream(std::istream& in, fastq_read& r);

    virtual void trim(int left, int right) override;
    virtual void write_to(std::ostream& out) const override;

    std::string const& get_qual() const;

private:

    std::string qual_;

};

struct sam_read
    : public read {

    sam_read();
    void init(std::string const& line);

    virtual void trim(int left, int right) override;
    virtual void trim(int left, int right, read& other, int o_left, int o_right) override;
    virtual void write_to(std::ostream& out) const override;

    bool paired() const;

private:

    int         flag_;
    std::string rname_;
    int         pos_;
    int         mapq_;
    std::string cigar_;
    std::string rnext_;
    int         pnext_;
    int         tlen_;
    std::string qual_;
    std::string optional_;

};

struct unexpected_nt : std::exception {};


int get_nt_number(char nt);
std::string reverse_complement(std::string const& seq);
bool seq_kmer(std::string const& seq, size_t k, size_t begin, u64& kmer);

} // namespace atlas
