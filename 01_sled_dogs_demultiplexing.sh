#Demultiplexing
#pool by pool
fastqpath=Pools/
  fastqfile=TOG-8IUJM-Fecal1_S1_L001

AdapterRemoval --file1 ${fastqpath}/${fastqfile}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal1 \
--barcode-list /Barcodes/Barcodes1.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile=TOG-8IUJM-Fecal2_S2_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal2 \
--barcode-list /Barcodes/Barcodes2.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile3=TOG-8IUJM-Fecal3_S3_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile3}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile3}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal3 \
--barcode-list /Barcodes/Barcodes3.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile4=TOG-8IUJM-Fecal4_S4_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile4}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile4}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal4 \
--barcode-list /Barcodes/Barcodes4.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile5=TOG-8IUJM-Fecal5_S5_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile5}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile5}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal5 \
--barcode-list /Barcodes/Barcodes5.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile6=TOG-8IUJM-Fecal6_S6_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile6}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile6}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal6 \
--barcode-list /Barcodes/Barcodes6.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile7=TOG-8IUJM-Fecal7_S7_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile7}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile7}_R2_001.fastq.gz \
--basename Demultiplex2/Fecal7 \
--barcode-list /Barcodes/Barcodes7.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile8=TOG-8IUJM-Fecal8_S8_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile8}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile8}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal8 \
--barcode-list /Barcodes/Barcodes8.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile9=TOG-8IUJM-Fecal9_S9_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile9}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile9}_R2_001.fastq.gz \
--basename /Demultiplex2/Fecal9 \
--barcode-list /Barcodes/Barcodes9.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile10=TOG-8IUJM-Saliva1_S10_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile10}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile10}_R2_001.fastq.gz \
--basename /Demultiplex2/Saliva1 \
--barcode-list /Barcodes/Barcodes10.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#Second run
fastqfile11=TOG-QIAJN-Saliva2_S1_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile11}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile11}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva2 \
--barcode-list /Barcodes/Barcodes11.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile12=TOG-QIAJN-Saliva3_S2_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile12}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile12}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva3 \
--barcode-list /Barcodes/Barcodes12.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile13=TOG-QIAJN-Saliva4_S5_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile13}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile13}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva4 \
--barcode-list /Barcodes/Barcodes13.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile14=TOG-QIAJN-Saliva5_S6_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile14}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile14}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva5 \
--barcode-list /Barcodes/Barcodes14.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile15=TOG-QIAJN-Saliva6_S7_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile15}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile15}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva6 \
--barcode-list /Barcodes/Barcodes15.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile16=TOG-QIAJN-Saliva7_S8_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile16}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile16}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva7 \
--barcode-list /Barcodes/Barcodes16.txt --barcode-mm-r1 2 --barcode-mm-r2 2

fastqfile17=TOG-QIAJN-Saliva8_S9_L001
AdapterRemoval --file1 ${fastqpath}/${fastqfile17}_R1_001.fastq.gz \
--file2 ${fastqpath}/${fastqfile17}_R2_001.fastq.gz \
--basename /Demultiplex3/Saliva8 \
--barcode-list /Barcodes/Barcodes17.txt --barcode-mm-r1 2 --barcode-mm-r2 2

#Remove unnecesary files
rm /Demultiplex2/*.discarded
rm /Demultiplex2/*.settings
rm /Demultiplex2/*.singleton.truncated
rm /Demultiplex2/*unidentified*
  
rm /Demultiplex3/*.discarded
rm /Demultiplex3/*.settings
rm /Demultiplex3/*.singleton.truncated
rm /Demultiplex3/*unidentified*
  
#Rename files
cd /Demultiplex2/
samples=$(ls *.pair1.truncated | sed 's/.pair1.truncated//')
for sample in ${samples}
do
mv ${sample}.pair1.truncated ${sample}.1.fq
done
mv ${sample}.pair2.truncated ${sample}.2.fq

#second run
cd /Demultiplex3/
  samples=$(ls *.pair1.truncated | sed 's/.pair1.truncated//')
for sample in ${samples}; do
mv ${sample}.pair1.truncated ${sample}.1.fq
mv ${sample}.pair2.truncated ${sample}.2.fq
done

#REMOVE PRIMERS
# Get samples name
ls *.1.fq | sed "s/\.1\.fq//g" > samples
#Pool1
cd /Pool1
ls *_1_trimmed.fq | sed "s/\_1\_trimmed\.fq//g" > samples
ls *_2_trimmed.fq | sed "s/\_2\_trimmed\.fq//g" > samples2

#Pool2
cd /Pool2
ls *_1_trimmed.fq | sed "s/\_1\_trimmed\.fq//g" > samples2

#Remove primers
#pool1
for sample in $(cat samples)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o /Trimmed1/${sample}_1_trimmed.fq -p /Trimmed1/${sample}_2_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1

echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCANCAG \
--discard-untrimmed \
-o /Trimmed2/${sample}_2rev_trimmed.fq -p /Trimmed2/${sample}_1rev_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1
done

#Pool2

for sample in $(cat samples)
do
echo "On sample: $sample"
cutadapt -e 0.15 -g ^CTANGGGNNGCANCAG -G ^GACTACNNGGGTATCTAAT \
--discard-untrimmed \
-o /Trimmed1/${sample}_1_trimmed.fq -p /Trimmed1/${sample}_2_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1

echo "On sample: $sample"
cutadapt -e 0.15 -g ^GGACTACNNGGGTATCTAAT -G ^CCTANGGGNNGCANCAG \
--discard-untrimmed \
-o /Trimmed2/${sample}_2rev_trimmed.fq -p /Trimmed2/${sample}_1rev_trimmed.fq \
${sample}.1.fq ${sample}.2.fq \
>> cutadapt_primer_trimming_stats.txt 2>&1
done