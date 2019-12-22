'''
make_fasta.py - make fastas of spliced genes + proteins in mtMinus

requires output from bedtools getfasta w/ fullHeader enabled
'''

from Bio import SeqIO
import argparse
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='make fastas of mtminus genes', 
        usage='python3.5 make_fasta.py [options]')

    parser.add_argument('-f', '--fasta', required=True,
                        type=str, help='Fasta from bedtools getfasta')
    parser.add_argument('-g', '--gff', required=True,
                        type=str, help='GFF file')
    parser.add_argument('-n', '--nuc', required=False,
                        action='store_true', help='return nuc sequences instead \
                        of protein')
    parser.add_argument('-o', '--outdir', required=True,
                        type=str, help='Dir to write fastas to')

    args = parser.parse_args()

    return args.fasta, args.gff, args.nuc, args.outdir

def get_coords(gff):
    with open(gff, 'r') as f:
        gff_all = [line.split() for line in f if 'gene' in line]
    fwd_coords = [l for l in gff_all if l[6] == '+']
    rev_coords = [l for l in gff_all if l[6] == '-']
    # format - god this is ugly
    fwd_coords = [[int(l[3])-1, int(l[4]), l[-1].split(';')[0].split('=')[1]] for
            l in fwd_coords]
    rev_coords = [[int(l[3])-1, int(l[4]), l[-1].split(';')[0].split('=')[1]] for
            l in rev_coords]
    # output looks like [0, 100, 'genename']
    return fwd_coords, rev_coords


def get_sequences(fasta):
    fwd_genes = [s for s in SeqIO.parse(fasta, 'fasta') if '(+)' in s.id]
    rev_genes = [s for s in SeqIO.parse(fasta, 'fasta') if '(-)' in s.id]
    fwd_genes_dict, rev_genes_dict = {}, {}
    for record in fwd_genes:
        fwd_genes_dict[str(record.id)] = record
    for record in rev_genes:
        rev_genes_dict[str(record.id)] = record
    return fwd_genes, rev_genes, fwd_genes_dict, rev_genes_dict

def get_cds_dicts(fwd_coords, rev_coords, fwd_genes, rev_genes):
    fwd_cds = {}
    rev_cds = {}
    print('generating fwd gene dict')
    for record in tqdm(fwd_genes):
        start, end = [int(val) for val in str(record.id).lstrip('mtMinus:').rstrip('(+)').split('-')]
        # oh boy this will be slow
        for gene_start, gene_end, gene_name in fwd_coords:
            if gene_start in [77918, 81014]: # the only two fwd one-CDS genes, lol
                fwd_cds[gene_start] = [record.id]
            elif start == gene_start and end == gene_end:
                continue
            elif start >= gene_start and end < gene_end:
                if gene_start not in fwd_cds:
                    fwd_cds[gene_start] = []
                fwd_cds[gene_start].append(record.id)
            elif start > gene_start and end <= gene_end:
                if gene_start not in fwd_cds:
                    fwd_cds[gene_start] = []
                fwd_cds[gene_start].append(record.id)

    print('generating rev gene dict')
    for record in tqdm(rev_genes[::-1]):
        start, end = [int(val) for val in str(record.id).lstrip('mtMinus:').rstrip('(-)').split('-')]
        for gene_start, gene_end, gene_name in rev_coords:
            if gene_start == 345009: # the only rev one-CDS gene
                rev_cds[gene_start] = [record.id]
            elif start == gene_start and end == gene_end:
                continue
            elif start >= gene_start and end < gene_end:
                if gene_start not in rev_cds:
                    rev_cds[gene_start] = []
                rev_cds[gene_start].append(record.id)
            elif start > gene_start and end <= gene_end:
                if gene_start not in rev_cds:
                    rev_cds[gene_start] = []
                rev_cds[gene_start].append(record.id)

    return fwd_cds, rev_cds

def write_proteins(fwd_coords, rev_coords, fwd_genes_dict, rev_genes_dict,
        fwd_cds, rev_cds, nuc, outdir):
    for start, end, gene_name in tqdm(fwd_coords):
        fname_out = outdir + gene_name + '.fasta'
        with open(fname_out, 'w') as f:
            f.write('>' + gene_name + ':{s}-{e}(+)'.format(s=start, e=end) + '\n')
            sequence_out = ''
            for cds in fwd_cds[start]:
                if nuc:
                    sequence_out += str(fwd_genes_dict[cds].seq)
                elif not nuc:
                    try:
                        assert len(fwd_genes_dict[cds].seq) % 3 == 0
                    except:
                        print(len(fwd_genes_dict[cds]))
                        print(fname_out)
                    sequence_out += str(fwd_genes_dict[cds].seq.translate())
            f.write(sequence_out + '\n')
    for start, end, gene_name in tqdm(rev_coords):
        fname_out = outdir + gene_name + '.fasta'
        with open(fname_out, 'w') as f:
            f.write('>' + gene_name + ':{s}-{e}(-)'.format(s=start, e=end) + '\n')
            sequence_out = ''
            for cds in rev_cds[start]:
                if nuc:
                    sequence_out += str(rev_genes_dict[cds].seq)
                elif not nuc:
                    try:
                        assert len(rev_genes_dict[cds].seq) % 3 == 0
                    except:
                        print(len(rev_genes_dict[cds]))
                        print(fname_out)
                    sequence_out += str(rev_genes_dict[cds].seq.translate())
            f.write(sequence_out + '\n')

def main():
    fasta, gff, nuc, outdir = args()
    if not outdir[-1] == '/':
        outdir += '/'
    fwd_coords, rev_coords = get_coords(gff)
    fwd_genes, rev_genes, fwd_genes_dict, rev_genes_dict = get_sequences(fasta)
    fwd_cds, rev_cds = get_cds_dicts(fwd_coords, rev_coords, fwd_genes, rev_genes)
    write_proteins(fwd_coords, rev_coords, fwd_genes_dict, rev_genes_dict,
            fwd_cds, rev_cds, nuc, outdir)

if __name__ == '__main__':
    main()

        

