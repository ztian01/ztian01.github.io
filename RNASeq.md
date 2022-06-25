# RNASeq Implementation
**Languages:** Shell, R <br>
**Softwares / packages:** sratoolkit, FastQC, MultiQC, Trimmomatic, fastp, HISAT2, Samtools, featureCounts, DESeq2 <br>
**Data:**  Himes, Blanca E., et al. "RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive
 gene that modulates cytokine function in airway smooth muscle cells." PloS one 9.6 (2014): e99625.
## Download SRA files
A SRR_Acc_list.txt file was downloaded from SRA Run Selector. Sratoolkit was used to download .sra files. <br>

```<language>
mkdir ../sraData
/home/ligw/projects/Bioinfo/Tools/sratoolkit.3.0.0-centos_linux64/bin/prefetch --option-file ../SRR_Acc_List.txt -p -O /home/ligw/projects/Bioinfo/learnRNASeq/20220623/sraData &
```

![sraDownload](RNASeq%20Implementation/20220623/Records/sraDownload.png)

## Fastp
Fastp