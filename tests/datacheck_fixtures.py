"""Published-layout fixtures shared by adapter and real Nextflow smoke tests."""
import csv
from pathlib import Path

from test_check_transcript_preservation import annotation, extend

AUDIT_HEADER = ['chromosome', 'strand', 'end', 'splice_site', 'gene_ids',
                'original_selected_endpoint', 'protected_endpoint', 'annotation_endpoint', 'outcome', 'blocked_by']
STATS_HEADER = ['gene_id', 'transcript_id', 'transcript_location', 'five_prime_extension',
                'three_prime_extension', 'five_prime_support', 'five_prime_read_count',
                'three_prime_support', 'three_prime_read_count']
FOLDERS = ('raw_annotation', 'extended_transcripts', 'extended_transcripts_stats',
           'fixed_extensions_before_partial_protection', 'fixed_extensions',
           'partial_protection', 'natural_overlaps')


def write_run(root, empty_negative=False):
    root = Path(root)
    for folder in FOLDERS:
        (root / folder).mkdir(parents=True, exist_ok=True)
    raw_parts = []
    for strand, label, tid, gene, offset in (
            ('+', 'positive', 't_plus.1', 'g_plus', 0), ('-', 'negative', 't_minus.1', 'g_minus', 400)):
        raw = annotation(identifier=tid, strand=strand).replace('g1', gene)
        # Keep negative exons in transcript order, as the published raw does.
        if strand == '-':
            lines = raw.splitlines(True)
            exons = [line for line in lines if '\texon\t' in line]
            rows = iter(reversed(exons))
            raw = ''.join(next(rows) if '\texon\t' in line else line for line in lines)
        final = extend(raw)
        final = final.replace(';transcript_type=protein_coding', ';transcript_type=protein_coding;leap_tag=LEAP_5,LEAP_3')
        def shifted(text):
            result = []
            for line in text.splitlines():
                f = line.split('\t')
                f[3], f[4] = str(int(f[3])+offset), str(int(f[4])+offset)
                result.append('\t'.join(f)+'\n')
            return ''.join(result)
        raw, final = shifted(raw), shifted(final)
        gene_row = lambda a, b: '1\ttest\tgene\t{}\t{}\t.\t{}\t.\tID={};gene_id={}\n'.format(a,b,strand,gene,gene)
        if not (empty_negative and strand == '-'):
            raw_parts.append(gene_row(100+offset,350+offset)+raw)
        final = '' if empty_negative and strand == '-' else gene_row(90+offset,360+offset)+final
        (root/'extended_transcripts'/('1_'+label+'_extended.gff')).write_text(final)
        with (root/'extended_transcripts_stats'/('1_'+label+'_stats.tsv')).open('w',newline='') as handle:
            writer = csv.DictWriter(handle,fieldnames=STATS_HEADER,delimiter='\t');writer.writeheader()
            if final:
                writer.writerow(dict(gene_id=gene,transcript_id=tid,transcript_location='1:{}-{}'.format(90+offset,360+offset),
                    five_prime_extension=10,three_prime_extension=10,five_prime_support='reads',five_prime_read_count=7,
                    three_prime_support='reads',three_prime_read_count=7))
        for kind in ('FE','LE'):
            left = (kind=='FE')==(strand=='+')
            a,b,baseline = (90+offset,150+offset,100+offset) if left else (300+offset,360+offset,350+offset)
            bed = '' if not final else '1\t{}\t{}\t{}\t0\t{}\t.\t.\t7\t{}\n'.format(a,b,baseline,strand,gene)
            (root/'fixed_extensions_before_partial_protection'/('test_1_'+label+'_'+kind+'_SCF.bed')).write_text(bed)
            (root/'fixed_extensions'/('1_'+label+'_'+kind+'_protected.bed')).write_text(bed)
        with (root/'partial_protection'/('1_'+label+'_partial_protection.tsv')).open('w',newline='') as handle:
            csv.writer(handle,delimiter='\t').writerow(AUDIT_HEADER)
        (root/'natural_overlaps'/('1_'+label+'_NO.txt')).write_text('')
    (root/'raw_annotation/current_annotation_1.gff3').write_text(''.join(raw_parts))
    (root/'reference.sizes').write_text('1\t2000\n')
    return root
