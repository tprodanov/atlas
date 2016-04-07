from amplicon_index import AmpliconIndex
import nt_string
import trim_read


def load_amplicons(f):
    return [(ampl.strip(), nt_string.reverse_complement(ampl.strip())) for ampl in f]


def create_amplicon_index(amplicons, k, f):
    ai = AmpliconIndex(k, amplicons)
    ai.save(f)


def main():
    # with open('../data/readsense/amplicons.txt') as f:
    #     amplicons = load_amplicons(f)
    with open('ai.txt', 'r') as f:
        ai = AmpliconIndex.load(f)
    k = ai.get_k()
    with open('../data/readsense/11_S17_L001_R1_001.fastq') as f:
        fastq = nt_string.read_fastq(f)

    guesses = [0] * 20
    for read in fastq:
        guesses[len(trim_read.possible_amplicons(read.seq, ai, k))] += 1
    with open('guesses.txt', 'w') as f:
        for g in guesses:
            f.write('%d\n' % g)


if __name__ == '__main__':
    main()
