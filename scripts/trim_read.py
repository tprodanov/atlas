import argparse
import nt_string


gap_penalty = -1.4
mismatch_penalty = -1.4
score_threshold = 80

h = [[-1e6] * 1000 for _ in range(1000)]


def trim_from_left_edge(read, amplicon):
    n = len(read)
    m = len(amplicon)
    # h = [[-1e6] * (n + 1) for _ in range(m + 1)]

    best_score = -1e6
    best_read_end = 0

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

    return best_score, 0, best_read_end


def trim_from_right_edge(read, amplicon):
    n = len(read)
    m = len(amplicon)
    # h = [[-1e6] * (n + 1) for _ in range(m + 1)]

    best_score = -1e6
    best_read_begin = n

    for i in range(m, -1, -1):
        h[i][n] = 0
        for j in range(n - 1, -1, -1):
            # without "if j"
            h[i][j] = h[i][j + 1] + gap_penalty
            if i < m:
                h[i][j] = max(h[i][j], h[i + 1][j] + gap_penalty)
                h[i][j] = max(h[i][j], h[i + 1][j + 1] +
                              (1.0 if amplicon[i] == read[j] else mismatch_penalty))
            if h[i][j] >= best_score:
                best_score = h[i][j]
                best_read_begin = j

    return best_score, best_read_begin, n


def trim_to_amplicon(read, amplicon):
    l_score, i1, j1 = trim_from_left_edge(read, amplicon)
    r_score, i2, j2 = trim_from_right_edge(read, amplicon)

    if l_score > r_score:
        return l_score, i1, j1
    else:
        return r_score, i2, j2


def possible_amplicons(read, amplicon_index, k):
    l = len(read)
    if l >= 2 * k:
        kmer_pos = [0,
                    k,
                    int(l / 2 - k / 2),
                    l - 2 * k,
                    l - k
        ]
    elif l >= k:
        kmer_pos = [0, l - k]
    else:
        return []

    all_amplicons = set()
    was_N = []
    for pos in kmer_pos:
        kmer = nt_string.count_kmer(read, pos, k)
        if kmer == -1:
            was_N.append(pos)
        else:
            for a in amplicon_index.get_all_amplicons(kmer):
                all_amplicons.add(a)
    if not all_amplicons and was_N:
        for pos in was_N:
            for kmer in nt_string.count_kmer_extend_n(read, pos, k):
                for a in amplicon_index.get_all_amplicons(kmer):
                    all_amplicons.add(a)
    return all_amplicons


def trim_reads1(reads, amplicons, amplicon_index, output, error_output):
    k = amplicon_index.get_k()
    for read in reads:
        best_score = -1
        best_i = 0
        best_j = 0
        read_seq = read.seq

        for ix in possible_amplicons(read_seq, amplicon_index, k):
            if ix > 0:
                score, i, j = trim_to_amplicon(read_seq, amplicons[ix - 1][0])
            else:
                score, i, j = trim_to_amplicon(read_seq, amplicons[-ix - 1][1])
            if score > best_score:
                best_score = score
                best_i, best_j = i, j

        if best_score > score_threshold:
            read.seq = read_seq[best_i:best_j]
            read.quality = read.quality[best_i:best_j]
            output.write(str(read))
        else:
            error_output.write(str(read))
