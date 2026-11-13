"""Small real alignment fixtures; no read names are treated as unique record IDs."""
from pathlib import Path
import subprocess


def alignments():
    def row(name, flag, chrom, start, cigar='5M3N5M', quality=60):
        return '\t'.join(map(str, (name, flag, chrom, start, quality, cigar, '*', 0, 0,
                                    'ACGTACGTAC', 'IIIIIIIIII', 'RG:Z:sample', 'NM:i:0')))
    return [row('shared', 0, '1', 100), row('reverse', 16, '1', 200),
            row('shared', 256, '1', 300), row('shared', 2064, '1', 400),
            row('duplicate', 1024, '2', 50), row('qc', 528, '2', 100, quality=0),
            row('mitochondrial', 0, 'MT', 20), row('unmapped', 4, '*', 0, '*', 0)]


def write_bam(path, rows=None, lengths=None, group='fixture'):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    lengths = lengths or {'1': 2000, '2': 2000, 'MT': 2000}
    sam = path.with_suffix('.sam')
    sam.write_text('@HD\tVN:1.6\tSO:unsorted\n' +
        ''.join('@SQ\tSN:{}\tLN:{}\n'.format(c, n) for c, n in lengths.items()) +
        '@RG\tID:sample\tSM:' + group + '\n@PG\tID:mapper\tPN:fixture\n' +
        ''.join(row + '\n' for row in (alignments() if rows is None else rows)))
    subprocess.run(['samtools', 'sort', '--no-PG', '-o', str(path), str(sam)], check=True, capture_output=True)
    index = Path(str(path) + '.csi')
    subprocess.run(['samtools', 'index', '-c', str(path), str(index)], check=True, capture_output=True)
    sam.unlink()
    return path, index


def sam_records(path, *flags):
    text = subprocess.check_output(['samtools', 'view', '--no-PG', *flags, str(path)], text=True)
    return text.splitlines()
