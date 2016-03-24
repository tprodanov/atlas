import argparse
import nt_string


gap_penalty = -1.4
mismatch_penalty = -1.4


def trim_from_left_edge(read, amplicon):
    n = len(read)
    m = len(amplicon)
    h = [[-1e6] * (n + 1) for _ in range(m + 1)]

    best_score = -1e6
    best_read_end = 1

    for i in range(m + 1):
        h[i][0] = 0
        for j in range(1, n + 1):
            # without "if j"
            h[i][j] = h[i][j - 1] + gap_penalty
            if i:
                h[i][j] = max(h[i][j], h[i - 1][j] + gap_penalty)
                h[i][j] = max(h[i][j], h[i - 1][j - 1] +
                              (1.0 if amplicon[i - 1] == read[j - 1] else mismatch_penalty))
            if h[i][j] >= best_score:
                best_score = h[i][j]
                best_read_end = j

    return best_score, read[:best_read_end]


def trim(read, amplicon, reversed_amplicon):
    score, trimmed_read = trim_from_left_edge(read, amplicon)

    reversed_read = ''.join(reversed(read))
    r_score, r_trimmed_read = trim_from_left_edge(reversed_read, reversed_amplicon)
    if r_score > score:
        score = r_score
        trimmed_read = reversed(r_trimmed_read)

    rc_read = nt_string.reverse_complement(read)
    rc_score, rc_trimmed_read = trim_from_left_edge(rc_read, amplicon)
    if rc_score > score:
        score = rc_score
        trimmed_read = nt_string.reverse_complement(rc_trimmed_read)

    c_read = ''.join(reversed(rc_read))
    c_score, c_trimmed_read = trim_from_left_edge(c_read, reversed_amplicon)
    if c_score > score:
        score = c_score
        trimmed_read = nt_string.reverse_complement(''.join(reversed(c_trimmed_read)))

    trimmed_read = ''.join(trimmed_read)
    return score, trimmed_read


def main():
    with open('test_in.fa') as fin:
        sequences = nt_string.read_fasta(fin)
    _, amplicon = sequences[0]
    _, read = sequences[1]
    r_amplicon = list(reversed(amplicon))

    score, trimmed = trim(read, amplicon, r_amplicon)
    print(read)
    print(trimmed)
    print(score)


if __name__ == '__main__':
    main()
