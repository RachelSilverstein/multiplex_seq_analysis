# Analysis software for arrayed multiplex variant sequencing


## Requirements:
* python 3.5.1
* conda 4.10.3
* bowtie 2.2.2
* trimgalore  0.4.1
* cutadapt  1.18
* FastQC v0.10.1
* samtools 1.4.1-1

### Special note for those using the Erisone cluster at MGH:
Load the required modules using:
```
module load python/3.5.1
module load bowtie2
module load trimgalore
module load samtools
module load cutadapt
```



## Installation:

`conda create -n "multiplex_seq_env" python=3.5.1`

`conda activate multiplex_seq_env`

`pip install -r requirements.txt`
