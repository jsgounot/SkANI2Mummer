## SkANI2Mummer

A simple python (>v.3 required) script to run SkANI first and then Mummer on pairs with distance >= `pd`. Can also run pairwise Mummer with SkANI prescreen. Minimal version of drep, without the clustering at the end.

### Install

```
mamba install -c conda-forge -c bioconda tqdm skani mummer4
```

### Usage

```
python skani2mummer.py --help
usage: python skani2mummer.py [-h] [-pd PD] [-s [S]] [-t T] [-q] flist outfile

positional arguments:
  flist       File containing fasta paths
  outfile     Outfile (compressed tsv.gz file)

options:
  -h, --help  show this help message and exit
  -pd PD      Primary distance / Min distance for Skani. 100 = Ignore SkANI and run all pairs (default: 90)
  -s [S]      Skani output file. Temporary file if not provided. (default: None)
  -t T        number of CPUs (default: 1)
  -q          Quiet TQDM (default: False)
  ```

You can run the script with example the files:

```
python skani2mummer.py example/flist.txt example/outfile.tsv.gz -t 16 -pd 85 -s example/skani.tsv
```