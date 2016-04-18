import argparse
import sys

from amplicon_index import AmpliconIndex
import nt_string
import trim_read


def load_amplicons(f):
    return [(amplicon.seq, nt_string.reverse_complement(amplicon.seq)) for amplicon in nt_string.read_fasta(f)]


def build_amplicon_index(args):
    amplicons = load_amplicons(args.amplicons)
    ai = AmpliconIndex(int(args.k), amplicons)
    ai.save(args.index)


def trim1(args):
    amplicons = load_amplicons(args.amplicons)
    ai = AmpliconIndex.load(args.index)
    reads = nt_string.read_fastq(args.fastq)
    trim_read.trim_reads1(reads, amplicons, ai, args.output, args.error)


def main():
    parser = argparse.ArgumentParser(prog='Amplicon Trimmer')
    subparsers = parser.add_subparsers(dest='command')

    parser_build = subparsers.add_parser('build', help='build amplicon index')
    parser_build.add_argument('-a', '--amplicons', type=argparse.FileType('r'), help='amplicons in fasta file')
    parser_build.add_argument('-i', '--index', type=argparse.FileType('w'), help='output: amplicon index')
    parser_build.add_argument('-k', type=int, default=13, help='kmer size')

    parser_trim1 = subparsers.add_parser('trim1', help='trim single reads')
    parser_trim1.add_argument('-a', '--amplicons', type=argparse.FileType('r'), help='amplicons in fasta file')
    parser_trim1.add_argument('-i', '--index', type=argparse.FileType('r'), help='amplicon index')
    parser_trim1.add_argument('-f', '--fastq', type=argparse.FileType('r'), help='fastq reads')
    parser_trim1.add_argument('-o', '--output',
                              type=argparse.FileType('w'), default=sys.stdout, help='output: trimmed reads')
    parser_trim1.add_argument('-e', '--error',
                              type=argparse.FileType('w'), default=sys.stderr, help='output: unmatched reads')

    # args = parser.parse_args('build -a ../data/11_S17/amplicons.fasta -i ../data/11_S17/ai13.txt'.split())
    args = parser.parse_args('trim1 -a ../data/11_S17/amplicons.fasta -i ../data/11_S17/ai13.txt'
                             ' -f ../data/11_S17/11_S17_1.fastq'
                             ' -o ../data/11_S17/11_S17_1_o.fastq'
                             ' -e ../data/11_S17/11_S17_1_e.fastq'.split())

    if args.command == 'build':
        build_amplicon_index(args)
    elif args.command == 'trim1':
        trim1(args)


if __name__ == '__main__':
    main()
