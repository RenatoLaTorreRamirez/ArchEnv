# ArchEnv
Analysis of ancient DNA from sedimentary samples
## Read preprocess
```
reads1="/path/to/reads1"
reads2="/path/to/reads2"
prefix=""
threads=30
mkdir -p ${prefix}_output

AdapterRemoval --file1 $reads1 --file2 $reads2 --mm 3 --minlength 30 --trimns --trimqualities --minquality 30 --basename ${prefix}_output/$prefix --threads $threads --collapse --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
### Merge collapsed and singleton reads
```
cat ${prefix}_output/${prefix}.collapsed ${prefix}_output/${prefix}.collapsed.truncated ${prefix}_output/${prefix}.singleton ${prefix}_output/${prefix}.singleton.truncated > ${prefix}_output/${prefix}.col.s.fq
```
### Clean intermediate files
```
rm ${prefix}_output/${prefix}.singleton
rm ${prefix}_output/${prefix}.singleton.truncated
rm ${prefix}_output/${prefix}.collapsed
rm ${prefix}_output/${prefix}.collapsed.truncated
rm ${prefix}_output/${prefix}.discarded
rm ${prefix}_output/${prefix}.settings
```
### Trim residual adapters and poly A from collapsed and singleton reads
```
fastq-grep --invert-match "AAAAA$" --mismatches="${prefix}_output/kmer_${prefix}.col.s.mis.fq" ${prefix}_output/${prefix}.col.s.fq > ${prefix}_output/kmer_${prefix}.col.s.fq
fastq-grep --invert-match "^TTTTT" --mismatches="${prefix}_output/kmer2_${prefix}.col.s.mis.fq" ${prefix}_output/kmer_${prefix}.col.s.fq > ${prefix}_output/kmer2_${prefix}.col.s.fq
seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ${prefix}_output/kmer2_${prefix}.col.s.fq \
| fastq-grep --invert-match "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" --mismatches="${prefix}_output/adap1_kmer2_${prefix}.col.s.mis.fq" > ${prefix}_output/adap1_kmer2_${prefix}.col.s.fq
seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^ACACTCTTTCCCTACACGACGCTCTTCCGATCT" ${prefix}_output/adap1_kmer2_${prefix}.col.s.fq \
| fastq-grep --invert-match "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" --mismatches="${prefix}_output/adap2_kmer2_${prefix}.col.s.mis.fq" > ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq
```
### Trim residual adapters and poly A from paired end reads
```
fastq-grep --invert-match "AAAAA$" --mismatches="${prefix}_output/kmer_${prefix}.pair1.truncated.mis" ${prefix}_output/${prefix}.pair1.truncated > ${prefix}_output/kmer_${prefix}.pair1.truncated
seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ${prefix}_output/kmer_${prefix}.pair1.truncated \
| fastq-grep --invert-match "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" --mismatches="${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.mis" > ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated
fastq-grep --invert-match "^TTTTT" --mismatches="${prefix}_output/kmer_${prefix}.pair2.truncated.mis" ${prefix}_output/${prefix}.pair2.truncated > ${prefix}_output/kmer_${prefix}.pair2.truncated
seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT$" ${prefix}_output/kmer_${prefix}.pair2.truncated \
| fastq-grep --invert-match "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" --mismatches="${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.mis" > ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated
```
Separate unpaired reads and add to collapsed/singleton file
```
fastq_pair ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated
if [ -s ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.single.fq ]; then
        cat ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.single.fq >> ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq;
fi
if [ -s ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.single.fq ]; then
        cat ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.single.fq >> ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq;
fi
```
### Trim and recover discarded reads (collapsed and singleton)
```
scraps=""
if [ -s ${prefix}_output/kmer_${prefix}.col.s.mis.fq ]; then
        seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AAAAA" ${prefix}_output/kmer_${prefix}.col.s.mis.fq \
        | fastq-grep --trim_after --trim_match "AAAAA" > ${prefix}_output/kmer_${prefix}.col.s.mis.trim.fq;
        scraps=1;
fi
if [ -s ${prefix}_output/kmer2_${prefix}.col.s.mis.fq ]; then
        seqtk seq -r ${prefix}_output/kmer2_${prefix}.col.s.mis.fq \
        | seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AAAAA" \
        | fastq-grep --trim_after --trim_match "AAAAA" \
        | seqtk seq -r > ${prefix}_output/kmer2_${prefix}.col.s.mis.trim.fq;
        scraps=1;
fi
if [ -s ${prefix}_output/adap1_kmer2_${prefix}.col.s.mis.fq ]; then
        fastq-grep --trim_after --trim_match "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ${prefix}_output/adap1_kmer2_${prefix}.col.s.mis.fq > ${prefix}_output/adap1_kmer2_${prefix}.col.s.mis.trim.fq;
        scraps=1;
fi
if [ -s ${prefix}_output/adap2_kmer2_${prefix}.col.s.mis.fq ]; then
        seqtk seq -r ${prefix}_output/adap2_kmer2_${prefix}.col.s.mis.fq \
        | seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^ACACTCTTTCCCTACACGACGCTCTTCCGATCT" \
        | fastq-grep --trim_after --trim_match "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" > ${prefix}_output/adap2_kmer2_${prefix}.col.s.mis.trim.fq;
        scraps=1;
fi

if [ "$scraps" = 1 ]; then
        cat ${prefix}_output/*col.s.mis.trim.fq >> ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq;
fi
```
### Trim and recover discarded reads (paired reads)
```
pscraps=""
if [ -s ${prefix}_output/kmer_${prefix}.pair1.truncated.mis ]; then
        seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AAAAA" ${prefix}_output/kmer_${prefix}.pair1.truncated.mis \
        | fastq-grep --trim_after --trim_match "AAAAA" > ${prefix}_output/kmer_${prefix}.pair1.truncated.mis.trim;
        pscraps=1;
fi
if [ -s ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.mis ]; then
        fastq-grep --trim_after --trim_match "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.mis > ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.mis.trim;
        pscraps=1;
fi
if [ -s ${prefix}_output/kmer_${prefix}.pair2.truncated.mis ]; then
        seqtk seq -r ${prefix}_output/kmer_${prefix}.pair2.truncated.mis \
        | seqkit grep --use-regexp --by-seq --invert-match --only-positive-strand --pattern "^AAAAA" \
        | fastq-grep --trim_after --trim_match "AAAAA" \
        | seqtk seq -r > ${prefix}_output/kmer_${prefix}.pair2.truncated.mis.trim;
        pscraps=1;
fi
if [ -s ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.mis ]; then
        seqtk seq -r ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.mis \
        | fastq-grep --trim_after --trim_match "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" \
        | seqtk seq -r > ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.mis.trim;
        pscraps=1;
fi

if [ "$pscraps" = 1 ]; then
        cat ${prefix}_output/*kmer*.truncated.mis.trim \
        | seqkit fx2tab \
        | sed 's/^/S_/' \
        | seqkit tab2fx >> ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq;
fi
```
### Remove low complexity sequences
```
sga preprocess --remove-adapter-fwd=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --remove-adapter-rev=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --dust-threshold=4 --min-length=30 ${prefix}_output/adap2_kmer2_${prefix}.col.s.fq --out=${prefix}_output/adap2_kmer2_${prefix}.col.s.pp.fq
bbmap/dedupe.sh threads=$threads in=${prefix}_output/adap2_kmer2_${prefix}.col.s.pp.fq out=${prefix}_output/${prefix}.unpaired.dedup.fq

sga preprocess --remove-adapter-fwd=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --remove-adapter-rev=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --dust-threshold=4 --min-length=30 --pe-mode=1 ${prefix}_output/adap1_kmer_${prefix}.pair1.truncated.paired.fq ${prefix}_output/adap1_kmer_${prefix}.pair2.truncated.paired.fq --out=${prefix}_output/adap1_kmer_${prefix}.pair12.truncated.paired.pp.fq
bbmap/dedupe.sh interleaved=true threads=$threads in=${prefix}_output/adap1_kmer_${prefix}.pair12.truncated.paired.pp.fq out=${prefix}_output/adap1_kmer_${prefix}.pair12.truncated.paired.pp.dedup.fq
bbmap/reformat.sh in=${prefix}_output/adap1_kmer_${prefix}.pair12.truncated.paired.pp.dedup.fq out1=${prefix}_output/${prefix}.pair1.dedup.fq out2=${prefix}_output/${prefix}.pair2.dedup.fq
```
### Clean intermediate files
```
rm ${prefix}_output/adap*
rm ${prefix}_output/kmer*
rm ${prefix}_output/*.truncated
rm ${prefix}_output/${prefix}.col.s.fq
```
## Kraken filter
### Initial microorganism classification using plusPF database
```
mkdir -p ${prefix}_output/01_${prefix}_KPF
kraken2 --db /path/to/KRAKEN2_PLUSPF --threads $threads --minimum-base-quality 20 --report ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.c.sequences --unclassified-out ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.unc.sequences --output ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.output ${prefix}_output/${prefix}.unpaired.dedup.fq 2>&1 \
| tee ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.log
kraken2 --paired --db /path/to/KRAKEN2_PLUSPF --threads $threads --minimum-base-quality 20 --report ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired#.c.sequences --unclassified-out ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired#.unc.sequences --output ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired.output ${prefix}_output/${prefix}.pair1.dedup.fq ${prefix}_output/${prefix}.pair2.dedup.fq 2>&1 \
| tee ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired.log
combine_kreports.py --reports ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.report ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired.report --output ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.combined.report --only-combined --no-headers
```
### Eukaryothe pathogen classification using EuPathDB database
```
mkdir -p ${prefix}_output/02_${prefix}_KEP
kraken2 --db /path/to/KRAKEN2_EuPath --threads $threads --minimum-base-quality 20 --report ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.c.sequences --unclassified-out ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.unc.sequences --output ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.output ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.unpaired.unc.sequences 2>&1 \
| tee ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.log
kraken2 --paired --db /path/to/KRAKEN2_EuPath --threads $threads --minimum-base-quality 20 --report ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired#.c.sequences --unclassified-out ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired#.unc.sequences --output ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired.output ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired_1.unc.sequences ${prefix}_output/01_${prefix}_KPF/${prefix}.plusPF.paired_2.unc.sequences 2>&1 \
| tee ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired.log
combine_kreports.py --reports ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.report ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired.report --output ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.combined.report --only-combined --no-headers
```
### Full NCBI NT classification using 2023 DB
```
mkdir -p ${prefix}_output/03_${prefix}_KFN
kraken2 --db /path/to/KRAKEN2_NCBI_NT --threads $threads --minimum-base-quality 20 --report ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.c.sequences --unclassified-out ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.unc.sequences --output ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.output ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.unpaired.unc.sequences 2>&1 \
| tee ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.log
kraken2 --paired --db /path/to/KRAKEN2_NCBI_NT --threads $threads --minimum-base-quality 20 --report ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired.report --confidence 0.05 --minimum-hit-groups 3 --report-minimizer-data --classified-out ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired#.c.sequences --unclassified-out ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired#.unc.sequences --output ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired.output ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired_1.unc.sequences ${prefix}_output/02_${prefix}_KEP/${prefix}.EuPath.paired_2.unc.sequences 2>&1 \
| tee ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired.log
combine_kreports.py --reports ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.unpaired.report ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.paired.report --output ${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.combined.report --only-combined --no-headers
```
## Reduce report
```
report=${prefix}_output/03_${prefix}_KFN/${prefix}.NCBI_NT.combined.report
total_lines=$(cat $report | wc -l)

extract_lvl () {
        tax_lvl=$1
        tax_name=$2
        total_lvl=$(cat $report | grep -P "\t$tax_lvl\t" | wc -l)
        lvl_lim=$(cat $report | grep -Pn "\t$tax_lvl\t" | grep -n . | grep "$tax_name" | cut -d':' -f1)
        lim1=$(cat $report | grep -Pn "\t$tax_lvl\t" | grep -n . | grep "$tax_name" | cut -d':' -f2)
        if [ "$lvl_lim" = "$total_lvl" ]; then
                cat $report \
                | sed -n "$lim1,${total_lines}p";
        else
                let lvl_lim2=lvl_lim+1
                lim2=$(cat $report | grep -Pn "\t$tax_lvl\t" | sed -n "${lvl_lim2}p" | cut -d':' -f1)
                let lim2=lim2-1
                cat $report \
                | sed -n "$lim1,${lim2}p"
        fi
}

if grep -q "R1" $report; then
        extract_lvl "R1" "cellular organisms" > ${report}.id1
        report=${report}.id1
fi
if grep -q "Eukaryota" $report; then
        extract_lvl "D" "Eukaryota" > ${report}.id2
        report=${report}.id2
fi
if grep -q "Viridiplantae" $report; then
        extract_lvl "K" "Viridiplantae" >> ${report}.euka
fi
if grep -q "Metazoa" $report; then
        extract_lvl "K" "Metazoa" >> ${report}.euka
fi
```

## Get combined summary report
```
project=""

mkdir -p Summary_${project}

cp *_output/03_*_KFN/*.euka Summary_${project}/
for i in $(ls Summary_${project}/*.euka | cut -d'/' -f2 | cut -d'.' -f1); do
        cat Summary_${project}/${i}.*.euka \
        | awk '{print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8}' > Summary_${project}/${i}.report;
done
files=$(ls Summary_${project}/*.report | tr '\n' ' ' | sed '$s/ $/\n/')
kraken-biom --max D --min SS --fmt tsv --output_fp Summary_${project}/Summary_${project}_D_S_output.tsv $files

sed -i '1d' Summary_${project}/Summary_${project}_D_S_output.tsv
sed '1d' Summary_${project}/Summary_${project}_D_S_output.tsv \
| cut -f1 \
| taxonkit lineage -n -r > Summary_${project}/Summary_${project}_D_S_output.lineages
taxonkit reformat --output-ambiguous-result --format "{K};{k};{p};{c};{o};{f};{g};{s};{t}" Summary_${project}/Summary_${project}_D_S_output.lineages \
| csvtk -H -t cut -f 1,4,3,5 \
| csvtk -H -t sep -f 4 -s ';' -R \
| csvtk add-header -t -n taxid,rank,name,superkingdom,kingdom,phylum,class,order,family,genus,species,subspecies \
| paste - Summary_${project}/Summary_${project}_D_S_output.tsv > Summary_${project}/Summary_${project}_D_S_full.tsv
```
