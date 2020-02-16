#!/usr/bin/env python

from random import choice
from dark.reads import DNARead

seq1 = ('NTTAAATCTGCTTGTTAGGCCTAAAAATAATGTCTTAC'
        'TTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTT')
seq2 = ('NATTTGAAGAGGACAAAGATGCTTGGAGCTAAGAAAAA'
        'AGATTCTGAAACACACTGAAACTTGTAATTATAAGTAC')

quality1 = 'Q' * len(seq1)
quality2 = 'Q' * len(seq2)

randomSeq = ''.join(choice('ACGT') for _ in range(200))

with open('query-1.fastq', 'w') as fp:
    print('@query\n%s\n+\n%s' % (seq1, quality1), file=fp)

with open('query-2.fastq', 'w') as fp:
    print('@query\n%s\n+\n%s' % (seq2, quality2), file=fp)

read2 = DNARead('query', seq2).reverseComplement()

with open('reference-matching.fasta', 'w') as fp:
    print('>reference-matching\n%s\n%s\n%s' %
          (seq1, randomSeq, read2.sequence), file=fp)

referenceNonMatching = ''.join(choice('ACGT') for _ in range(
    len(seq1) + len(randomSeq) + len(seq2)))

with open('reference-non-matching.fasta', 'w') as fp:
    print('>reference-non-matching\n%s' % referenceNonMatching, file=fp)


PAIRED = 0x1
UNMAP = 0x4
MUNMAP = 0x8
PAIR1 = 0x40
PAIR2 = 0x80

with open('query.sam', 'w') as fp:
    print('@HD\tVN:1.0\tSO:unsorted', file=fp)
    print('@SQ\tSN:ref-0\tLN:150', file=fp)
    print('\t'.join(map(str, (
        'query', PAIRED | UNMAP | MUNMAP | PAIR1, '*', 0, 0, '*', '*', 0, 0,
        seq1, quality1))), file=fp)
    print('\t'.join(map(str, (
        'query', PAIRED | UNMAP | MUNMAP | PAIR2, '*', 0, 0, '*', '*', 0, 0,
        seq2, quality2))), file=fp)
