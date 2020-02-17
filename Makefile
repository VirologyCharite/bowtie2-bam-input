all: data match-fastq match-sam

data:
	./make-data.py
	samtools view -O BAM query.sam > query.bam
	bowtie2-build --quiet reference-matching.fasta reference-matching
	bowtie2-build --quiet reference-non-matching.fasta reference-non-matching

match-fastq:
	bowtie2 --quiet --no-head -x reference-matching -1 query-1.fastq -2 query-2.fastq

non-match-fastq:
	bowtie2 --quiet -x reference-non-matching -1 query-1.fastq -2 query-2.fastq > non-match.sam
	cat non-match.sam
	samtools view -O BAM non-match.sam > non-match.bam
	bowtie2 --quiet --no-head -x reference-matching -b non-match.bam --align-paired-reads

match-sam:
	bowtie2 --quiet --no-head -x reference-matching -b query.sam --align-paired-reads

match-bam:
	bowtie2 --quiet --no-head -x reference-matching -b query.bam --align-paired-reads

clean:
	rm -f query-1.fastq query-2.fastq query.bam query.sam non-match.bam non-match.sam
	rm -f reference-matching.* reference-non-matching.*
	rm -f *~
