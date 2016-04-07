import math
import nt_string
import pickle


class AmpliconIndex:

    def __init__(self, k, amplicons):
        self.__k = k
        self.__modulo = 4 ** (k - 1)

        self.__kmer_dict = dict()
        amplicon_number = 0
        for amplicon, r_amplicon in amplicons:
            self.__add_amplicon(amplicon, amplicon_number)
            self.__add_amplicon(r_amplicon, amplicon_number)
            amplicon_number += 1

    def __add_amplicon(self, amplicon, amplicon_number):
        amplicon_iter = iter(amplicon)
        i = 1
        kmer = 1
        for nt in amplicon_iter:
            kmer = 4 * kmer + nt_string.nt_number[nt]
            i += 1
            if i == self.__k:
                break

        for nt in amplicon_iter:
            kmer = (kmer % self.__modulo) * 4 + nt_string.nt_number[nt]

            if kmer in self.__kmer_dict:
                self.__kmer_dict[kmer].append(amplicon_number)
            else:
                self.__kmer_dict[kmer] = set(amplicon_number)

    def get_all_amplicons(self, kmer):
        return self.__kmer_dict[kmer] if kmer in self.__kmer_dict else []

    def get_k(self):
        return self.__k

    def save(self, f):
        f.write('%d\n' % self.__k)
        for kmer, entries in self.__kmer_dict.items():
            f.write('%d:' % kmer)
            f.write('%s\n' % ' '.join(str(entry) for entry in entries))

    @staticmethod
    def load(f):
        self = AmpliconIndex.__new__(AmpliconIndex)
        self.__k = int(next(f))
        self.__modulo = 4 ** (self.__k - 1)
        self.__kmer_dict = dict()

        for line in f:
            kmer, entries = line.split(':')
            self.__kmer_dict[int(kmer)] = [int(entry) for entry in entries.split(' ')]
        return self

    def print(self):
        print(self.__kmer_dict)
