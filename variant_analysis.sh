# Downloading the reference genome
 wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz

# Indexing the reference genome with bwa
bwa index GCA_000195955.2_ASM19595v2_genomic.fna

# Obtaining the SRA file
prefetch SRR29748661

# Converting the SRA file into fastq files
 fastq-dump --split-files SRR29748661

# Performing QC with fastqc
mkdir -p qc_results
mv SRR29748661_1.fastq SRR29748661_2.fastq ~/qc_results
fastqc -o qc_results *.fastq
xdg-open SRR29748661_1_fastqc.html SRR29748661_2_fastqc.html
sudo apt-get install trim-galore
trim_galore --paired  SRR29748661_1.fastq SRR29748661_2.fastq

# Aligning the fastq reads to the reference genome with bwa
 bwa mem GCA_000195955.2_ASM19595v2_genomic.fna   SRR29748661_1_val_1.fq  SRR29748661_2_val_2.fq > aligned_reads.sam

# Converting the sam file to a bam file with samtools
samtools view -Sb aligned_reads.sam > aligned_reads.bam

# Sorting the bam file with samtools
samtools sort aligned_reads.bam -o aligned_reads_sorted.bam

# Indexing the sorted bam file with samtools
samtools index aligned_reads_sorted.bam

# Calling the variants using bcftools
bcftools mpileup -Ou -f GCA_000195955.2_ASM19595v2_genomic.fna aligned_reads_sorted.bam | bcftools call -mv -Ov -o SRR29748661_variants.vcf
less -S SRR29748661_variants.vcf

# Filtering the variants
bcftools filter -s LowQual -e '%QUAL<20 || DP<10' SRR29748661_variants.vcf -Ov -o filtered_variants.vcf
less -S filtered_variants.vcf

# Downloading the database for M. tuberculosis and annotating
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
java -jar snpEff/snpEff.jar databases | grep 'Mycobacterium_tuberculosis'
java -jar snpEff/snpEff.jar download Mycobacterium_tuberculosis_h37rv
java -jar snpEff/snpEff.jar Mycobacterium_tuberculosis_h37rv filtered_variants.vcf > annotated_variants.vcf

# Summarizing the annotation data
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' annotated_variants.vcf | \
> cut -d '|' -f 1,2,3,4,5 > summary_variants.txt
nano summary_variants.txt

# Analysis
sudo apt install vcftools
bcftools stats filtered_variants.vcf

# I perform Visualisation with IGV

# Functional anotation with snpEff to understand the biological impact of each variant
snpeff ann -csvStats snpeff_stats.csv -csvStats snpeff_annotated.vcf > snpeff_summary.txt

# Filtering variants
bcftools filter -i 'INFO/DP > 20' filtered_variants.vcf -o high_depth_variants.vcf # Filtering by depth

bcftools filter -i 'INFO/AC > 1' filtered_variants.vcf -o high_ac_variants.vcf # Filtering by allele count of variants

cut -d'|' -f3 summary_variants.txt | sort | uniq # To get the different allele impacts
