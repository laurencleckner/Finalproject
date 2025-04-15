# Final project

---
Our goal is to compare gut samples of healthy patients and breast cancer patients. Here, we can quantify the differences in bacteria present, potenitally identifying differencies in quantity and species.

We were planning on looking at microbiota from breast cancer tissue samples as well, but we were unable to do that becuase the samples were single and not paired. Due to this, we will just be looking at the gut microbiota data. 

All of our data should be paired and amplicon. 

For this project we will be using:

QUIIME2 --> to quantify the amount of bacteria present in each sample. This is available in the HPC.

---
## **Step 1: Downloading the data**
1. Close the repositroy to get the accession_list.txt file
```
git clone repositorylink
```
2. Create a new directory for the project and move the accession-list.txt file into the directory
```
mkdir -p microbiome_project/fasta
mv accession_list.txt microbiome_project/
```
3. Go into the new directory and check if the file is there
```
cd microbiome_project
ls
```
4. Now we will need to download the FASTQ files from SRA. To do this we will need to create a file and save the script and run bash.

```
vi prep_data.sh
```
```
#!/bin/bash

# Create folders & clear any previous versions
mkdir -p fastq
echo "sample-id,absolute-filepath,direction" > manifest.csv
echo -e "sample-id\tPatientStatus" > metadata.tsv

# Loop over your accession list
while read line; do
    accession=$(echo $line | cut -f1 -d' ')
    status=$(echo $line | cut -f2 -d' ')

    echo "Downloading $accession..."
    fastq-dump --split-files $accession --outdir fastq

    filepath=$(realpath fastq/${accession}_1.fastq)  # single-end assumed
    echo "$accession,$filepath,forward" >> manifest.csv
    echo -e "$accession\t$status" >> metadata.tsv
done < accession_list.txt
```
5. Load the SRA tool kit
```
module load sra-toolkit/2.10.9
```
6. Automate the download of all your raw sequencing data from NCBI using the list of accession numbers. This is result in 2 files: manifest.csv and metadata..tsv
```
bash prep_data.sh
```
You can now check these files to make sure that they look correct
```
head manifest.csv
cat metadata.tsv
```
In the metadta file we will need to change some names of samples
```
sed -i 's/Breast/Cancer/g' metadata.tsv
```
---
## **Step 2: Data Quality Control**

1. Run fastqc
```
module load FastQC
mkdir fastqc_raw_data_out
fastqc -t 4 fastq/*.fastq -o fastqc_raw_data_out
```
2. Run multiqc (will not work if you activate conda environment first, in such case `conda deactivate` before running multiqc)
```
module load python/3.8.6
multiqc --dirs fastqc_raw_data_out --filename multiqc_raw_data.html
```
3. Add, commit and push file to github repository. You can also transfer files using GLOBUS connect.
```
git add "multiqc_raw_data.html"
git commit -m "Adding multiqc results"
git push origin main
```
4. Go to web browser github and download `multiqc_raw_data.html`
5. Open file in web browser
---
## **Step 3: Import Data into QIIME 2**
### Import the data:
1. We need to change our raw files names so that they can be read by qiime
- Create directory for renamed files
```
mkdir -p casava_reads
```
  - Open the vi editor
```
vi rename.sh
```
  - Type `I` to edit and write the following:
```
#!/bin/bash
for file in fastq/*_1.fastq; do
    base=$(basename "$file" _1.fastq)
    
    cp "fastq/${base}_1.fastq" "casava_reads/${base}_S1_L001_R1_001.fastq"
    cp "fastq/${base}_2.fastq" "casava_reads/${base}_S1_L001_R2_001.fastq"
done
```
- Run rename.sh script using:
```
bash rename.sh
```
#### OOPS some sequences are missing (do not do this step if you were not in Jorges office)
- Download missing sequences.
```
vi tissue_accessions
```
- Copy and paste missing accessions
```
vi download_sra.sh
```
```
#!/bin/bash
# Use prefetch to download all files
module load sra-toolkit
while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR 
   fastq-dump --gzip $SRR --split-files -O fastq_files/
done < tissue_accessions.txt
```
- We realized the tissue sequences were single not paired, we decided to remove them from the analysis
2. Let's create the following directory
```         
mkdir reads_qza
```
3. Now we will create the environment
```
module load anaconda3        
conda activate qiime2-amplicon-2024.2
```
4. Compress fastq files --> we needed to do this becasue QUIIME2 would not accept them when unzipped
```
gzip SRR*
```
When you do this step make sure that you cd into the casava_reads directory 

5. Feed qiime the raw reads (make sure to be in the microbiome_proect directory when doing this step)
```         
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path casava_reads \
  --output-path reads_qza/reads.qza \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt
```
- cassava format, files follow a pattern `SampleID_L001_R1_001.fastq.gz`
6. remove primer sequences from reads, these are the primers used to enrich for a specific locus, e.g.:16S, COI, etc
```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads_qza/reads.qza \
  --p-cores 4 \
  --p-front-f ACGCGHNRAACCTTACC \
  --p-front-r ACGGGCRGTGWGTRCAA \
  --p-discard-untrimmed \
  --p-no-indels \
  --o-trimmed-sequences reads_qza/reads_trimmed.qza
```
New code, relaxed:

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads_qza/reads.qza \
  --p-cores 4 \
  --p-front-f ACGCGHNRAACCTTACC \
  --p-front-r ACGGGCRGTGWGTRCAA \
  --p-match-adapter-wildcards \
  --p-error-rate 0.1 \
  --o-trimmed-sequences reads_qza/reads_trimmed.qza
```

7. Visualize your data now
```         
qiime demux summarize \
  --i-data reads_qza/reads_trimmed.qza \
  --o-visualization reads_qza/reads_trimmed_summary.qzv
```
- Add, commit and push `reads_qza/reads_trimmed_summary.qzv` to git
```
git add "reads_qza/reads_trimmed_summary.qzv"
git commit -m "Adding reads summary results"
git push origin main
```
- Download and open `reads_trimmed_summary.qzv` in https://view.qiime2.org/
---

## **Step 4: Denoising with DADA2**
#### This is importnat to filter out low-quality sequences, chimeras, and artifacts. This allows for more accurate data and downstream analysis reliability.
1. Ask for some computer resources
```
salloc --mem=100G --time=6:00:00 --cpus-per-task=32 
```
### Important issues here, cancer samples have reads 150 bp long, whereas healthy samples are 300 bp

-We are running dada separately.
```
 awk 'NR>2 && $3==0 && $4==0
 && $5==0 && $6==0 && $7==0 && $8==0 && $9==0 {print $1}' metadata.tsv > accessions_150.txt

awk 'NR>2 && !($3==0 && $4=
=0 && $5==0 && $6==0 && $7==0 && $8==0 && $9==0) {print $1}' metadata.tsv > accessions_300.txt
```
- Now convert to qiime format
```
echo -e "#SampleID" > samples_150.tsv
cat accessions_150.txt >> samples_150.tsv

echo -e "#SampleID" > samples_300.tsv
cat accessions_300.txt >> samples_300.tsv
```
- Split reads
```
qiime demux filter-samples \
  --i-demux reads.qza \
  --m-metadata-file samples_150.tsv \
  --o-filtered-demux reads_150.qza

qiime demux filter-samples \
  --i-demux reads.qza \
  --m-metadata-file samples_300.tsv \
  --o-filtered-demux reads_300.qza
```

2. Run DADA2 to denoise and generate ASVs:

```
vi dada.slurm
```
- Type I and paste
```
#!/bin/bash
#SBATCH --job-name=dada2_denoise
#SBATCH --output=dada2_denoise.out
#SBATCH --error=dada2_denoise.err
#SBATCH --time=06:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=drandal1@svsu.edu

# Load conda and activate the QIIME2 environment
module load anaconda3
conda activate qiime2-amplicon-2024.2

# Run DADA2 denoising on cancer reads first
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads_150.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 140 \
  --p-trunc-len-r 140 \
  --p-n-threads 32 \
  --o-table table_150.qza \
  --o-representative-sequences repseqs_150.qza \
  --o-denoising-stats stats_150.qza

# Now run DADA2 denoising on healthy reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads_300.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 220 \
  --p-n-threads 32 \
  --o-table table_300.qza \
  --o-representative-sequences repseqs_300.qza \
  --o-denoising-stats stats_300.qza

```
- Exit : esc :wq
- Run the script

```
sbatch dada.slurm
```


2. Summarize the feature table:
```bash
qiime feature-table summarize \
  --i-table table_150.qza \
  --o-visualization table_150.qzv

qiime feature-table summarize \
  --i-table table_300.qza \
  --o-visualization table_300.qzv

qiime metadata tabulate \
  --m-input-file stats_150.qza \
  --o-visualization stats_150.qzv

qiime metadata tabulate \
  --m-input-file stats_300.qza \
  --o-visualization stats_300.qzv
```


3. Merge groups
```
qiime feature-table merge \
  --i-tables table_150.qza \
  --i-tables table_300.qza \
  --o-merged-table table_merged.qza

qiime feature-table merge-seqs \
  --i-data repseqs_150.qza \
  --i-data repseqs_300.qza \
  --o-merged-data repseqs_merged.qza
```
## **Step 5: Taxonomic Assignment**
1. Assign taxonomy using a pre-trained classifier (e.g., SILVA):
```bash
#Download database to compare to
wget https://data.qiime2.org/2023.9/common/silva-138-99-nb-classifier.qza

qiime feature-classifier classify-sklearn \
   --i-classifier silva-138-99-nb-classifier.qza \
   --i-reads repseqs_merged.qza \
   --p-n-jobs 8 \
   --output-dir taxa
```
---
### Filtering resultant table
1. Filter out rare ASVs
```         
qiime feature-table filter-features \
  --i-table table_merged.qza \
  --p-min-frequency 2 \
  --p-min-samples 1 \
  --o-filtered-table dada_table_filt.qza
```
- replace 2 and 1 with appropiate frequencies and observe the changes.
2. Filter out contaminant and unclassified ASVs
```         
qiime taxa filter-table \
  --i-table dada_table_filt.qza \
  --i-taxonomy taxa/classification.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table dada_table_filt_contam.qza
```
3. Subset and summarize filtered table
```         
qiime feature-table summarize \
  --i-table dada_table_filt_contam.qza \
  --o-visualization dada_table_filt_contam_summary.qzv
```
- Git add, commit and push `dada_table_filt_contam_summary.qzv` Download and open in https://view.qiime2.org/

4. Copy final table to current directory
```         
cp dada_table_filt_contam.qza dada_table_final.qza
```
5. Produce final table summary
```         
qiime feature-table filter-seqs \
  --i-data repseqs_merged.qza \
  --i-table dada_table_final.qza  \
  --o-filtered-data rep_seqs_final.qza
```
```         
qiime feature-table summarize \
  --i-table dada_table_final.qza \
  --o-visualization dada_table_final_summary.qzv
```

- Git add, commit and push `dada_table_final_summary.qzv` Download and open in https://view.qiime2.org/

## **Visualization**
2. Visualize taxonomic composition:
```bash
qiime taxa barplot \
    --i-table table_merged.qza \
    --i-taxonomy taxa/classification.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv
```
