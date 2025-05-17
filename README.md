# Project-1-Variant-Calling-Pipeline-FASTQ-VCF-
🧬 Variant Calling Pipeline: FASTQ → Annotated VCF  This project demonstrates a complete variant calling pipeline starting from raw sequencing reads in FASTQ format to annotated variants in VCF using open-source bioinformatics tools.

## 📌 Goal
To perform variant discovery from a human Whole Exome/Genome Sequencing dataset using standard best-practice tools, from alignment to annotation, using the **GRCh38** human reference genome.

## 🛠️ Tools & Software

| Tool      | Purpose                              |
|-----------|--------------------------------------|
| FastQC    | Quality check for FASTQ files        |
| MultiQC   | Aggregate QC reports                 |
| BWA       | Read alignment to reference          |
| SAMtools  | Convert, sort, and index BAM files   |
| GATK      | Reference indexing (dictionary)      |
| BCFtools  | Variant calling and VCF manipulation|
| SnpEff    | Variant annotation                   |

## 🗂️ Project Folder Structure

variant-calling-project/
├── data/
│   ├── fastq/           # Raw reads (FASTQ)
│   ├── ref/             # Reference genome and indexes
│   └── known/           # Known variants (optional)
├── results/
│   ├── align/           # SAM/BAM files
│   ├── variants/        # VCF and BCF files
│   └── annotation/      # SnpEff outputs
├── scripts/             # Custom scripts
└── logs/                # Log files (optional). 
snpEFF folder was outside the variant-calling-project/

## ⚙️ Setup & Installation
### Install core tools:
sudo apt update  
sudo apt install fastqc multiqc bwa samtools bcftools openjdk-11-jre-headless -y

### Install GATK:
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip  
unzip gatk-4.5.0.0.zip  
sudo mv gatk-4.5.0.0 /opt/gatk  
echo 'export PATH=$PATH:/opt/gatk' >> ~/.bashrc  
source ~/.bashrc  

### Install SnpEff:
cd ~  
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip  
unzip snpEff_latest_core.zip

# 🧪 Step-by-Step Pipeline
## 1️⃣ Download FASTQ Files
cd data/fastq  
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016004/SRR016004_1.fastq.gz  
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR016/SRR016004/SRR016004_2.fastq.gz  

## 2️⃣ Download and Prepare Reference Genome
cd data/ref  
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  

### Indexing
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa  
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa  
gatk CreateSequenceDictionary -R Homo_sapiens.GRCh38.dna.primary_assembly.fa  

Note: If gatk throws a Python error, make sure to modify the shebang to use python3.

## 3️⃣ Align Reads to Reference (BWA)
cd results/align  
bwa mem -t 2 ../../data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa ../../data/fastq/SRR016004_1.fastq.gz ../../data/fastq/SRR016004_2.fastq.gz > aln-pe.sam  

## 4️⃣ Convert SAM to BAM, Sort, and Index
samtools view -bS aln-pe.sam > aln-pe.bam  
samtools sort aln-pe.bam -o aln-pe.sorted.bam  
samtools index aln-pe.sorted.bam  

## 5️⃣ Variant Calling (bcftools)
cd ../variants  
### Create raw BCF
bcftools mpileup -f ../../data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa ../align/aln-pe.sorted.bam -Ou -o raw.bcf
### Call variants
bcftools call -mv -Ob -o variants.bcf raw.bcf  
### Convert to VCF
bcftools view variants.bcf > variants.vcf  

## 6️⃣ Variant Annotation with SnpEff
### Download GRCh38.86 Database
cd ~/snpEff  
java -jar snpEff.jar download GRCh38.86  

### Annotate Variants
cd ~/variant-calling-project/results/annotation  
java -Xmx4g -jar ../../../snpEff/snpEff.jar GRCh38.86 ../variants/variants.vcf > annotated.vcf    

## 📝Generate Summary Report
java -Xmx4g -jar ../../../snpEff/snpEff.jar GRCh38.86 ../variants/variants.vcf -stats stats.html -csvStats stats.csv > annotated.vcf  



