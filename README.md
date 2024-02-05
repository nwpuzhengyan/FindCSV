# FindCSV
## Introduction

Structural variations (SVs) play a significant role in genetic diseases and evolutionary mechanisms. Extensive research has been conducted over the past decade to detect simple structural variations, leading to the development of well-established detection methods. However, recent studies have highlighted the potentially greater impact of complex structural variations (CSVs) on individuals compared to simple structural variations. Despite this, the field still lacks precise detection methods specifically designed for complex structural variations. Therefore, the development of a highly efficient and accurate detection method is of utmost importance. In response to this need, we propose a novel method called FindCSV, which leverages deep learning techniques and consensus sequences to enhance the detection of SVs using long-read sequencing data. FindCSV demonstrates superior performance in detecting both complex and simple structural variations, as supported by experimental data.

---
## Installation
```
git clone https://github.com/nwpuzhengyan/FindCSV.git
```
---
## Dependence
    1. python3
	2. pysam
	3. cigar
	4. numpy
	5. pyfaidx
	6. copy
	7. time
	8. argparse
	9. PIL
	10. pytorch
	11. torchvision
	12. os
	13. swalign
	
---
## Running
The sorted bam files from NGMLR, Minimap and Minimap2 are all be used as input sorted bam. The input reference.fa and reference.fa of bam file must be the same one.
```
cd dist
./FindCSV <input sorted bam> <input reference.fa>	
```
---
## Convert SV regions into images

FindCSV can convert SV regions into images and store the images in this folder SV_into_image. If users want to check the converted images, they can look for them in folder SV_into_image. The name of each image consists of the chromosome name and SV start and end position.

---
## Output format
The output format is as follows. CHROM is chromosome name. POS is the SV start position. ID is the SV name. REF is the reference sequence and ALT is the alternate sequence. QUAL is the quality of SV and FILTER means filter status. INFO is the basic information of SV.
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10780	FindCSV.INS.1	G	GAACACATGCTAGCGCGTCCGGGGGTGGAGGCGATAGCGCAGGCGCAGAGAGCGCCGCGCC	.	PASS	SVTYPE=INS;SVLEN=61;END=10780;BP=chr1:10780-10780|;
chr1	30893	FindCSV.DEL.1	catttctctctatctcatttctctctctctcgctatct	c	.	PASS	SVTYPE=DEL;SVLEN=-37;END=30930;BP=chr1:30839-30930|;
```
---
## Contact
For advising, bug reporting and requiring help, please contact yan.zheng@nwpu-bioinformatics.com.
