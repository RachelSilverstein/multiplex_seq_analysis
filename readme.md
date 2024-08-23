# Analysis software for arrayed multiplex variant sequencing

## Requirements:
* python 3.x.x
* bowtie 2.2.2
* trimgalore  0.4.1
* cutadapt  1.18
* FastQC v0.10.1
* samtools 1.4.1-1


## Installation:

```
$ python -m venv multiplex_seq_env
$ source multiplex_seq_env/bin/activate
$ pip install -r requirements.txt
```


### Special note for those using the Erisone/two cluster at MGH:
Load the required modules using the following. For Eristwo SKIP loading python.
```
module load python/3.5.1
module load bowtie2
module load trimgalore
module load samtools
module load cutadapt
```

Sometimes loading the above modules will switch the version of python. If you get an
error try activating the venv again and the normal verion of python should be
reloaded.
You can check whether the right version of python is running in your environment
using `python --version`. Currently the correct version is the default python 3.8.8
and the incorrect version is python 3.6.6.
