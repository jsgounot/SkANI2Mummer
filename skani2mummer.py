import glob
import argparse
import os
import csv
import gzip
import logging
import tempfile
import subprocess
import shutil
import concurrent.futures
import itertools

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

try:
    from tqdm import tqdm
    itqdm = True
except ImportError:
    itqdm = False

def run_cmdline(cmdline):
    try:
        r = subprocess.run(cmdline, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        print(f"{err} {err.stderr.decode('utf8')}")
        exit(1)

def check_fasta(path, clean_path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rb') as f_in:
            with open(clean_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return clean_path
    return path

def get_dnadiff_prc(report_value):
    left_index = report_value.index('(') + 1
    return float(report_value[left_index:-2])

def parse_dnadiff_report(fname):
    subres = {}

    with open(fname) as f:
        for line in f:
            if line.startswith('AlignedBases'):
                line = line.strip().split()
                prc1, prc2 = get_dnadiff_prc(line[1]), get_dnadiff_prc(line[2])
                subres['mummer_prc_aligned1'], subres['mummer_prc_aligned2'] = prc1, prc2
                
            elif line.startswith('AvgIdentity'):
                line = line.strip().split()
                prc1, prc2 = float(line[1]), float(line[2])
                subres['mummer_avg_identity1'], subres['mummer_avg_identity2'] = prc1, prc2

    return subres

def run_mummer_pair(line, outdir=None):
    s1, s2 = line['Ref_file'], line['Query_file']

    # Run delta first
    tmp_dir = tempfile.TemporaryDirectory()
    s1 = check_fasta(s1, os.path.join(tmp_dir.name, 's1.fasta'))
    s2 = check_fasta(s2, os.path.join(tmp_dir.name, 's2.fasta'))

    #logger.info(f'Run mummer for {pair}')
    prefix = os.path.join(tmp_dir.name, 'prefix')
    cmdline = f'nucmer -p {prefix} {s1} {s2}'
    run_cmdline(cmdline)
    
    #logger.info(f'Run dnadiff for {pair}')
    delta = prefix + '.delta'
    prefix = os.path.join(tmp_dir.name, 'prefix.dnadiff')
    cmdline = f'dnadiff -d {delta} -p {prefix}'
    run_cmdline(cmdline)

    report = prefix + '.report'
    rdata = parse_dnadiff_report(report)

    if outdir:
        idx = line['idx']
        outfile = os.path.join(outdir, str(idx) + '.tar.gz')
        transform = f'--transform="s|{tmp_dir.name[1:]}|./{idx}|"'
        cmdline = f'tar --exclude="*.fasta" {transform} -czvf {outfile} {tmp_dir.name}'
        run_cmdline(cmdline)

    line.update(rdata)
    return line

def list_from_skani_other(args):
    ntf = tempfile.NamedTemporaryFile()
    outfile = args.s if args.s else ntf.name

    cmdline = f'skani dist --ql {args.flist} --rl {args.a} -o {outfile} -t {args.t}'
    run_cmdline(cmdline)

    dist = float(args.d)

    with open(outfile) as f:
        header = next(f).strip().split('\t')
        reader = csv.DictReader(f, delimiter='\t', fieldnames=header)
        lines = [{'idx': idx, ** line} for idx, line in enumerate(reader) if float(line['ANI']) >= dist]

    header.insert(0, 'idx')
    return header, lines

def list_from_skani_self(args):
    ntf = tempfile.NamedTemporaryFile()
    outfile = args.s if args.s else ntf.name

    cmdline = f'skani triangle -l {args.flist} -o {outfile} -t {args.t} --sparse'
    run_cmdline(cmdline)

    dist = float(args.d)

    with open(outfile) as f:
        header = next(f).strip().split('\t')
        reader = csv.DictReader(f, delimiter='\t', fieldnames=header)
        lines = [{'idx': idx, ** line} for idx, line in enumerate(reader) if float(line['ANI']) >= dist]

    header.insert(0, 'idx')

    return header, lines

def list_from_flist_other(args):
    with open(args.flist) as f:
        f1 = [line.strip() for line in f]

    with open(args.a) as f:
        f2 = [line.strip() for line in f]

    header = ['Ref_file', 'Query_file']
    return header, [
        {'idx': idx, 'Ref_file': p1, 'Query_file': p2}
        for idx, (p1, p2) in enumerate(itertools.product(f1, f2))
        ]

def list_from_flist_self(args):
    with open(args.flist) as f:
        fnames = [line.strip() for line in f]

    header = ['idx', 'Ref_file', 'Query_file']
    return header, [
        {'idx': idx, 'Ref_file': p1, 'Query_file': p2}
        for idx, (p1, p2) in enumerate(itertools.combinations(fnames, 2))
        ]

def _read_done_pairs(outfile):
    done = set()
    with gzip.open(outfile, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            done.add((row['Ref_file'], row['Query_file']))
    return done

def main(args):
    dist, a = float(args.d), args.a
    if dist == 100:
        header, lines = list_from_flist_other(args) if a else list_from_flist_self(args)
    elif 0 < dist < 100:
        header, lines = list_from_skani_other(args) if a else list_from_skani_self(args)
    else:
        raise Exception(f'd value must be 0 < d <= 100')

    full_header = header + ['mummer_prc_aligned1', 'mummer_avg_identity1',
                            'mummer_prc_aligned2', 'mummer_avg_identity2']
    checkpoint = int(args.c)

    resuming = False
    if checkpoint > 0 and os.path.exists(args.outfile):
        done = _read_done_pairs(args.outfile)
        if done:
            n_before = len(lines)
            lines = [l for l in lines if (l['Ref_file'], l['Query_file']) not in done]
            logger.info('Resuming: %d/%d pairs already done, %d remaining', len(done), n_before, len(lines))
            resuming = True

    logger.info('%d mummer pairs to run', len(lines))

    threads = int(args.t)
    use_tqdm = itqdm and not args.q
    fun = lambda line: run_mummer_pair(line, outdir=args.m)
    if args.m and not os.path.isdir(args.m): os.makedirs(args.m)

    n_done = 0
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(fun, line) for line in lines]
        iterator = concurrent.futures.as_completed(futures)
        if use_tqdm:
            iterator = tqdm(iterator, total=len(lines))

        if checkpoint > 0:
            out_mode = 'at' if resuming else 'wt'
            with gzip.open(args.outfile, out_mode) as f:
                writer = csv.DictWriter(f, full_header, delimiter='\t')
                if not resuming:
                    writer.writeheader()
                buffer = []
                for future in iterator:
                    buffer.append(future.result())
                    n_done += 1
                    if not use_tqdm and n_done % 1000 == 0:
                        logger.info('%d/%d pairs done', n_done, len(lines))
                    if len(buffer) >= checkpoint:
                        writer.writerows(buffer)
                        f.flush()
                        buffer = []
                if buffer:
                    writer.writerows(buffer)
        else:
            for future in iterator:
                results.append(future.result())
                n_done += 1
                if not use_tqdm and n_done % 1000 == 0:
                    logger.info('%d/%d pairs done', n_done, len(lines))

    if checkpoint == 0:
        with gzip.open(args.outfile, 'wt') as f:
            writer = csv.DictWriter(f, full_header, delimiter='\t')
            writer.writeheader()
            writer.writerows(results)

parser = argparse.ArgumentParser(
    prog='python skani2mummer.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('flist', help='File containing fasta paths')
parser.add_argument('outfile', help='Outfile (compressed tsv.gz file)')
parser.add_argument('-d', default=90, help='Min distance for Skani. 100 = Ignore SkANI and run all pairs')
parser.add_argument('-a', nargs='?', help='Target genomes files list. Otherwise do a self-comparison.')
parser.add_argument('-s', nargs='?', help='Skani output file. Temporary file if not provided.')
parser.add_argument('-m', nargs='?', help='Mummer output dir. All mummer files will be compressed and store in the directory.')
parser.add_argument('-t', default=1, help='number of CPUs')
parser.add_argument('-c', default=0, type=int, help='Save results to output every N pairs for crash recovery (0 = disabled)')
parser.add_argument('-q', action='store_true', help='Quiet TQDM')

args = parser.parse_args()
main(args)