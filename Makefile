all: data match-fastq match-sam

data:
	./make-data.py
	bowtie2-build --quiet reference-matching.fasta reference-matching
	bowtie2-build --quiet reference-non-matching.fasta reference-non-matching

match-fastq:
	bowtie2 --quiet --no-head -x reference-matching -1 query-1.fastq -2 query-2.fastq

non-match-fastq:
	bowtie2 --quiet --no-head -x reference-non-matching -1 query-1.fastq -2 query-2.fastq

match-sam:
	bowtie2 --quiet --no-head -x reference-matching -b query.sam
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --ff
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --fr
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --rf
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --align-paired-reads
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --align-paired-reads --ff
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --align-paired-reads --fr
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --align-paired-reads --rf

clean:
	rm -f query-1.fastq query-2.fastq query.sam
	rm -f reference-matching.* reference-non-matching.*
	rm -f *~
