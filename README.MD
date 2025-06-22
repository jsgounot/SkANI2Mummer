## SkANI2Mummer

A python (>v.3 required) script to run Mummer4 only for pairs with distances >= `d` from SkANI results. Can also run all pairs without SkANI screening by setting `d` = 100. Can be run against itself (default) or against another files list (`-a` option). 

All mummer files can be compressed and saved on the fly into a subdirectory with `-m` option. If set, all archives base name corresponds to the `idx` column in the output file.

### Install

```
mamba install -c conda-forge -c bioconda tqdm skani mummer4 python>=3.7
```

### Usage

```
python skani2mummer.py --help
usage: python skani2mummer.py [-h] [-d D] [-a [A]] [-s [S]] [-m [M]] [-t T] [-q] flist outfile

positional arguments:
  flist       File containing fasta paths
  outfile     Outfile (compressed tsv.gz file)

options:
  -h, --help  show this help message and exit
  -d D        Min distance for Skani. 100 = Ignore SkANI and run all pairs (default: 90)
  -a [A]      Target genomes files list. Otherwise do a self-comparison. (default: None)
  -s [S]      Skani output file. Temporary file if not provided. (default: None)
  -m [M]      Mummer output dir. All mummer files will be compressed and store in the directory. (default: None)
  -t T        number of CPUs (default: 1)
  -q          Quiet TQDM (default: False)
  ```

You can run the script with the example files:

```
python skani2mummer.py example/flist.txt example/outfile.tsv.gz -t 16 -pd 85 -s example/skani.tsv
```