The [bowtie2 docs](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
have a
[usage](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#usage)
section that says you can pass the queries as a BAM file, using

    -b <bam>

    Reads are unaligned BAM records sorted by read name. The
    --align-paired-reads and --preserve-tags options affect the way Bowtie
    2 processes records.

I (Terry) was trying get that to work. Turns out the file *must* be BAM
(not SAM) and you have to use `--align-paired-reads`.

## Workings

This repo has a Python script to generate some data and a `Makefile` to run
it.

```sh
$ make data
./make-data.py
bowtie2-build --quiet reference-matching.fasta reference-matching
bowtie2-build --quiet reference-non-matching.fasta reference-non-matching
```

Then we can do a match using `-1` and `-2` (some newlines inserted into the
output for readability):

```sh
$ make match-fastq
bowtie2 --quiet --no-head -x reference-matching -1 query-1.fastq -2 query-2.fastq
query	99	reference-matching	1	42	76M	=	277	352
NTTAAATCTGCTTGTTAGGCCTAAAAATAATGTCTTACTTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTT
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
AS:i:-1	XN:i:1	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:0N75	YS:i:-1	YT:Z:CP
query	147	reference-matching	277	42	76M	=	1	-352
GTACTTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTTTCTTAGCTCCAAGCATCTTTGTCCTCTTCAAATN	
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
AS:i:-1	XN:i:1	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:75N0	YS:i:-1	YT:Z:CP
```

The same queries can be in a BAM file. Try 

```$
$ make non-match-fastq
# No output from Bowtie2 because the queries do not match.
bowtie2 --quiet -x reference-non-matching -1 query-1.fastq -2 query-2.fastq > non-match.sam

# The FASTQ is in the SAM output, with flags indicating the lack of a match.
# Newlines inserted here for readability.
cat non-match.sam
@HD	VN:1.0	SO:unsorted
@SQ	SN:reference-non-matching	LN:352
@PG	ID:bowtie2	PN:bowtie2	VN:2.3.5.1	
    CL:"/home/linuxbrew/.linuxbrew/bin/../Cellar/bowtie2/2.3.5.1/bin/bowtie2-align-s
    --wrapper basic-0 --quiet -x reference-non-matching -1 query-1.fastq -2 query-2.fastq"
query	77	*	0	0	*	*	0	0	
NTTAAATCTGCTTGTTAGGCCTAAAAATAATGTCTTACTTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTT
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ	YT:Z:UP
query	141	*	0	0	*	*	0	0
NATTTGAAGAGGACAAAGATGCTTGGAGCTAAGAAAAAAGATTCTGAAACACACTGAAACTTGTAATTATAAGTAC
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ	YT:Z:UP

# Convert it to BAM.
samtools view -O BAM non-match.sam > non-match.bam

# Try again, using -b and the matching reference.
bowtie2 --quiet --no-head -x reference-matching -b non-match.bam --align-paired-reads
query	99	reference-matching	1	42	76M	=	277	352	
NTTAAATCTGCTTGTTAGGCCTAAAAATAATGTCTTACTTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTT
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
AS:i:-1	XN:i:1	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:0N75	YS:i:-1	YT:Z:CP
query	147	reference-matching	277	42	76M	=	1	-352
GTACTTATAATTACAAGTTTCAGTGTGTTTCAGAATCTTTTTTCTTAGCTCCAAGCATCTTTGTCCTCTTCAAATN
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
AS:i:-1	XN:i:1	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:75N0	YS:i:-1	YT:Z:CP
```

All good.

```sh
$ make clean
rm -f query-1.fastq query-2.fastq query.bam query.sam non-match.bam non-match.sam
rm -f reference-matching.* reference-non-matching.*
rm -f *~
```


