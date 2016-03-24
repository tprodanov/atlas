

def read_fasta(fin):
    name = next(fin).strip()
    current_sequence = []
    sequences = []
    for line in fin:
        if line[0] == '>':
            sequences.append((name, ''.join(current_sequence)))
            current_sequence = []
            name = line.strip()
        else:
            current_sequence.append(line.strip())
    sequences.append((name, ''.join(current_sequence)))
    return sequences


class UnexpectedLetter(Exception):
    pass


def complement_nt(nt):
    if nt == 'A':
        return 'T'
    if nt == 'C':
        return 'G'
    if nt == 'G':
        return 'C'
    if nt == 'T':
        return 'A'
    if nt == 'N':
        return 'N'
    raise UnexpectedLetter


def reverse_complement(seq):
    return ''.join(complement_nt(nt) for nt in reversed(seq))
