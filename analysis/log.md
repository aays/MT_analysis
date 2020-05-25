
## 8/12/2019

things to do:

- create a fasta of all the reinhardtii MT genes 
    - incl. FUS1, for completeness
- create fasta of reinhardtii MT gene proteins? (if doing blastp)
- include Volvox? talk this over with Rory
- figure out means of scripting reciprocal best blast

need to check for/get
- gene presence (Hamaji-like table)
    - should be done with tblastn
- gene synteny 
    - (ie need to retain positions on matched contigs as well, for downstream comparisons)
- give gene IDs if genes have been predicted
    - blastp against protein fasta file AND tblastn against genome
    - check if genomic coordinates match (check against gff) and get ID of predicted gene

## 21/12/2019

so this slipped me by with this hellish end-of-semester...

first up - making a FASTA of all the mtMinus genes

will then use the protein sequences of the predicted genes to tblastn
against the MT genes

doing this with `bedtools getfasta`

```bash
mkdir data/references
cd data/references

# get GFF
grep 'mtMinus' /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/annotation/concatenated_GFF/final.strict.GFF3 > mtMinus.gff

# get FASTA
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/mtMinus_ref.chromosome_6_and_mtMinus.fasta .

# make FASTAs
bedtools getfasta \
-fi data/references/mtMinus_ref.chromosome_5_and_mtMinus.fasta \
-bed data/references/mtMinus.gff \
-fullHeader -s \
-fo data/references/mtMinus_genes.fasta
```

an attempt at a script to create a translated 
version of all the fwd strands:

```bash
mkdir data/fastas
mkdir data/fastas-nuc

time python3.5 analysis/make_fasta.py \
--fasta data/references/mtMinus_genes.fasta \
--gff data/references/mtMinus.gff \
--nuc \
--outdir data/fastas-nuc
```

so this works - but some CDSs are not multiples of 3? - maybe just stick to nucleotide sequences

## 22/12/2019

combining all the nuc fastas into one:

```bash
touch data/fastas-nuc/fastas_all.fa
cat data/fastas-nuc/*.fasta >> data/fastas-nuc/fastas_all.fa;
```

tblastn with incerta genes as query and reinhardtii genes as subject -

(this seems to require a non-bgzipped file)

```bash
gunzip -c -d rory-data/Chlamydomonas_incerta.braker2.protein.fa.gz > rory-data/Chlamydomonas_incerta.braker2.protein.fa

time tblastn \
-query rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-subject data/fastas-nuc/fastas_all.fa \
-outfmt 6 \
-out data/test.out
```

naive run returns 50k matches, most (48k) of which have an evalue > 0.1 (49k hits > 0.01) - ran in 14 min

will repeat this with evalue thresh of 0.01 - sorting by bit score to find best matches

```bash
mkdir data/tblastn

time tblastn \
-query rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-subject data/fastas-nuc/fastas_all.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/tblastn/C_incerta_tblastn.tsv
```

things to do re: homology:
- repeat for other species
- figure out a way to do reciprocal hits
    - looking at incerta, several queries have 10+ hits, while one subject has 711 (!)
    - can do blastx with CDSs as query and incerta proteins as subject

repeating for other species:

```bash
for spec in Chlamydomonas_schloesseri Edaphochlamys_debaryana; do
    gunzip -c -d rory-data/${spec}.braker2.protein.fa.gz > rory-data/${spec}.braker2.protein.fa;
    gunzip -c -d rory-data/${spec}.braker2.CDS.fa.gz > rory-data/${spec}.braker2.CDS.fa;
    gunzip -c -d rory-data/${spec}.nuclear.fa.gz > rory-data/${spec}.nuclear.fa;
done
```

tblastn'ing:

```bash
time tblastn \
-query rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-subject data/fastas-nuc/fastas_all.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/tblastn/C_schloesseri_tblastn.tsv

time tblastn \
-query rory-data/Edaphochlamys_debaryana.braker2.protein.fa \
-subject data/fastas-nuc/fastas_all.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/tblastn/E_debaryana_tblastn.tsv
```

blastx'ing for reciprocal hits:

```bash
mkdir -p data/blastx

time blastx \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/blastx/C_incerta_blastx_rev.tsv

time blastx \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Edaphochlamys_debaryana.braker2.protein.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/blastx/E_debaryana_blastx_rev.tsv

time blastx \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/blastx/C_schloesseri_blastx_rev.tsv
```

going to look at these in R but first - next up - synteny

need to blastn full reinhardtii genes against species contigs to create presence/absence table
like in Hamaji 2018

to test for changes in synteny - create output file w/ matched contigs *and* where
on the contigs the matches are found - this is already mostly baked into `outfmt 6`

presence/absence:

```bash
mkdir -p data/scaffold-mapping

# first w/o evalue thresh to see what we get
time blastn \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Chlamydomonas_incerta.nuclear.fa \
-outfmt 6 \
-out data/scaffold-mapping/C_incerta_scaffold_hits.tsv
```

looks like e-value filtering not explicitly needed here - all matches have evals well below 0.01 -
will keep it in the command for consistency

```bash
time blastn \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Chlamydomonas_schloesseri.nuclear.fa \
-outfmt 6 -evalue 0.01 \
-out data/scaffold-mapping/C_schloesseri_scaffold_hits.tsv

time blastn \
-query data/fastas-nuc/fastas_all.fa \
-subject rory-data/Edaphochlamys_debaryana.nuclear.fa \
-outfmt 6 -evalue 0.01 \
-out data/scaffold-mapping/E_debaryana_scaffold_hits.tsv
```

now going to port these off the server for some RStudio fiddling + (manually?) making
the presence/absence table

tomorrow - sort out reciprocal best blast and reach out to rory

## 23/12/2019

porting reciprocal best blast outputs *and* gffs to local computer

should probably create filtered gffs only containing contigs with hits to reduce
the amt of space they take up

first, unzipping:

```bash
for spec in Chlamydomonas_incerta Chlamydomonas_schloesseri Edaphochlamys_debaryana; do
    gunzip -c -d rory-data/${spec}.braker2.gff3.gz > rory-data/${spec}.braker2.gff3;
done
```

going to just keep gene records, actually -

```bash
grep 'gene' rory-data/Chlamydomonas_incerta.braker2.gff3 > data/references/Chlamydomonas_incerta.genes.gff3
grep 'gene' rory-data/Chlamydomonas_schloesseri.braker2.gff3 > data/references/Chlamydomonas_schloesseri.genes.gff3
grep 'gene' rory-data/Edaphochlamys_debaryana.braker2.gff3 > data/references/Edaphochlamys_debaryana.genes.gff3
```

also - if checking reinhardtii genes against contigs - use tblastx instead of blastn! will
translate both query + subject contigs in all six frames

running tblastx for genome hits - 

```bash
time tblastx \
-query data/fastas-nuc/fastas-all.fa \
-subject rory-data/Chlamydomonas_incerta.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-out data/scaffold-mapping/C_incerta_scaffold_hits.tsv
```

took 21 min - queueing up for other two species

## 25/12/2019

do we see the same broken CDS problem in the MT+ annotation?
(which should be better than the MT- ones on paper)

```bash
grep 'chromosome_6' /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/annotation/concatenated_GFF/final.strict.GFF3 > chr6.gff
```

the coordinates for MT+ on chr6 were 298298-826737 (294645-826768 
technically it seems, given that those values overlap genes otherwise)

some quick awk filtering

```bash
awk -F '\t' '{ if(($4 >= 294645 && $5 <= 826768)) { print }}' data/references/chr6.gff > data/references/mtPlus.gff

getting only genes and CDS:

grep 'gene\|CDS' data/references/mtPlus.gff > data/references/mtPlus_genes.gff
```

## 28/12/2019

edit: forget the above - let's use Rob's already-extracted-and-named CDSs for this

getting Rob's fastas and the name translation dict:

```bash
cp -v /scratch/research/projects/chlamydomonas/mt_locus_recombination/analysis/cds-popgen/find_shared_genes/*CDS.fasta data/fastas-nuc

cp -v /scratch/research/projects/chlamydomonas/mt_locus_recombination/analysis/cds-popgen/find_shared_genes/*NameTranslation.txt data/
```

two data prep things to do here:
- convert fasta headers to gene names using NameTranslation dicts
- translate fastas from `fastas-nuc` into protein seqs in `fastas`

missing tab between column names in `NameTranslation` - needed a quick fix in python - same
goes for mtLimited name translation

code was messy but looks something like this:

```python
>>> with open('fastas-nuc/mtMinus_CDS_named.fasta', 'w') as f:
  2     for record in d_minus:
  3         outname = record.id
  4         for item in name_dict:
  5             if not 'alternate_ID' in item.keys():
  6                 if record.id == item['ch6_Ness_ID']:
  7                     outname = item['Common_Name']
  8                     break
  9                 elif record.id == item['mtMinus_Ness_ID']:
 10                     outname = item['Common_Name']
 11                     break
 12             else:
 13                 if record.id == item['alternate_ID']:
 14                     outname = item['Common_Name']
 15         if outname == record.id:
 16             print('wtf', record.id)
 17         f.write('>' + outname + '\n')
 18         f.write(str(record.seq) + '\n')
```

in short - replaced gene name with common name where there was one - otherwise,
retained alternate ID

now for translation:

```python
>>> from Bio import SeqIO
>>> from tqdm import tqdm
>>> d_plus = [s for s in SeqIO.parse('fastas-nuc/mtPlus_CDS_named.fasta', 'fasta')]
>>> d_minus = [s for s in SeqIO.parse('fastas-nuc/mtMinus_CDS_named.fasta', 'fasta')]
>>> with open('fastas/mtPlus_proteins.fasta', 'w') as f:
  2     for record in tqdm(d_plus):
  3         assert len(str(record.seq)) % 3 == 0
  4         f.write('>' + record.id + '\n')
  5         translated = str(record.seq.translate())
  6         f.write(translated + '\n')
100%|████████████████████████████████████████████████████████████████████████████| 114/114 [00:00<00:00, 2042.50it/s]
>>> with open('fastas/mtMinus_proteins.fasta', 'w') as f:
  2     for record in tqdm(d_minus):
  3         assert len(str(record.seq)) % 3 == 0
  4         f.write('>' + record.id + '\n')
  5         translated = str(record.seq.translate())
  6         f.write(translated + '\n')
100%|██████████████████████████████████████████████████████████████████████████████| 41/41 [00:00<00:00, 2258.02it/s]
>>> [r for r in d_plus if not str(r.seq.translate())[-1] == '*']
[]
>>> [r for r in d_minus if not str(r.seq.translate())[-1] == '*']
[]
```

seems all translated fine, and have stop codons at the end as expected - nice! 

back to blasting:
1. blastp - reinhardtii proteins vs other species proteins
    - reciprocal for all three species
2. tblastn - reinhardtii proteins vs other species genomes

blastp for orthology:

```bash
mkdir -p data/blastp

# 'forward' - species query, reinhardtii subject
time blastp \
-query rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-subject data/fastas/mtMinus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_incerta_fwd.tsv

time blastp \
-query rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-subject data/fastas/mtMinus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_schloesseri_fwd.tsv

time blastp \
-query rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-subject data/fastas/mtMinus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_schloesseri_fwd.tsv
# these three took 6 min each

# 'reverse' - reinhardtii query, species subject
time blastp \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_incerta_rev.tsv

time blastp \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_schloesseri_rev.tsv

time blastp \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/blastp/C_schloesseri_rev.tsv

```

after that, tblastn for gene content/synteny:

```bash
mkdir -p data/scaffold-mapping

time tblastn \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Chlamydomonas_incerta.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/scaffold-mapping/C_incerta_scaffold_hits.tsv

time tblastn \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/scaffold-mapping/C_schloesseri_scaffold_hits.tsv

time tblastn \
-query data/fastas/mtMinus_proteins.fasta \
-subject rory-data/Edaphochlamys_debaryana.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/scaffold-mapping/E_debaryana_scaffold_hits.tsv
# each takes about a minute to a minute 15
```

blasting w/ the plus proteins just to be sure that only shared genes show up:

```bash
time blastp \
-query rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-subject data/fastas/mtPlus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_incerta_fwd.tsv

time blastp \
-query rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-subject data/fastas/mtPlus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_schloesseri_fwd.tsv

time blastp \
-query rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-subject data/fastas/mtPlus_proteins.fasta \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_schloesseri_fwd.tsv

# 'reverse' - reinhardtii query, species subject
time blastp \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Chlamydomonas_incerta.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_incerta_rev.tsv

time blastp \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_schloesseri_rev.tsv

time blastp \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
-outfmt 6 \
-evalue 1e-10 \
-num_alignments 1 \
-max_hsps 1 \
-out data/plus-test/C_schloesseri_rev.tsv

# tblastn against the species genomes
time tblastn \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Chlamydomonas_incerta.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/plus-test/C_incerta_scaffold_hits.tsv

time tblastn \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Chlamydomonas_schloesseri.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/plus-test/C_schloesseri_scaffold_hits.tsv

time tblastn \
-query data/fastas/mtPlus_proteins.fasta \
-subject rory-data/Edaphochlamys_debaryana.nuclear.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/plus-test/E_debaryana_scaffold_hits.tsv
```

## 1/1/2020

next up - blasting predicted genes in the three unicells against the chlamy genome

primarily testing for translocations + genes that are not in chlamy but are in any or
all of these species

will be using protein sequences + tblastn again - though should probably
stick to genes that are on the contigs with hits in the reciprocal blast + tblastn

for incerta -
- C0033 (main)
- C0003 (just one gene - HDH1)

for schloesseri -
- C0045 (main)
- C0001 (97782)
- C0009 (PDK1)
- C0051 (SPP1C)

for debaryana -
- C0116 (most common)
- C0043 (second most common)
- C0001
- C0100
- C0077
- C0081
- C0012
- C0087
- C0125
- C0045

```bash
grep -e 'C00[03]3' rory-data/Chlamydomonas_incerta.braker2.gff3 > rory-data/Chlamydomonas_incerta.braker2.gff3.filtered
grep -e 'gene\|CDS' rory-data/Chlamydomonas_incerta.braker2.gff3.filtered | sponge rory-data/Chlamydomonas_incerta.braker2.gff3.filtered

awk -F '\t' '{ if($1 ~ /C0045|C0001|C0009|C0051/ ) { print } }' rory-data/Chlamydomonas_schloesseri.braker2.gff3
grep -e 'gene\|CDS' rory-data/Chlamydomonas_schloesseri.braker2.gff3.filtered | sponge rory-data/Chlamydomonas_schloesseri.braker2.gff3.filtered
```

debaryana has too many contigs for me to do this in awk over a crappy server connection - python it is

```python
>>> fname = 'rory-data/Edaphochlamys_debaryana.braker2.gff3'
>>> contigs = ['C0116', 'C0043', 'C0001', 'C0100', 'C0077', 'C0081', 'C0012', 'C0087', 'C0125', 'C0045']
>>> from tqdm import tqdm
  2 with open(fname + '.filtered', 'w') as f_out:
  3     with open(fname, 'r') as f_in:
  4         for line in tqdm(f_in):
  5             if not line.startswith('##'):
  6                 sp = line.split('\t')
  7                 contig, feature = sp[0], sp[2]
  8                 if contig in contigs and feature in ['gene', 'CDS']:
  9                     f_out.write(line)
815210it [00:01, 508896.67it/s]
```

starting with incerta, since it's probably simplest

need to make a filtered protein fasta with just proteins from these contigs

calling it a day for now, but read in contigs + gene names,
iterate through proteins w/ SeqIO, and if protein name in gff then print to new fasta

could probably write an actual generalized script for this and then run across all 
three species instead of coding this into the console every time

## 2/1/2019

attempting a python script to create these filtered fastas:

```bash
time python3.5 analysis/make_filtered_fasta.py \
--fasta rory-data/Chlamydomonas_incerta.braker2.protein.fa \
--gff rory-data/Chlamydomonas_incerta.braker2.gff3.filtered \
--outfile data/fastas/C_incerta.contigs.protein.fa
```

looks good - filtered down to 854 genes from 16957

repeating for other two species

```bash
time python3.5 analysis/make_filtered_fasta.py \
--fasta rory-data/Chlamydomonas_schloesseri.braker2.protein.fa \
--gff rory-data/Chlamydomonas_schloesseri.braker2.gff3.filtered \
--outfile data/fastas/C_schloesseri.contigs.protein.fa
# 1553 genes left, from 16268

time python3.5 analysis/make_filtered_fasta.py \
--fasta rory-data/Edaphochlamys_debaryana.braker2.protein.fa \
--gff rory-data/Edaphochlamys_debaryana.braker2.gff3.filtered \
--outfile data/fastas/E_debaryana.contigs.protein.fa
# 1804 genes left, from 20450
```

next up - tblastn these proteins against the chlamy genome

```
time tblastn \
-query data/fastas/C_incerta.contigs.protein.fa \
-subject data/references/chlamy.5.3.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/reinhardtii-tblastn/C_incerta_hits.tsv

time tblastn \
-query data/fastas/C_schloesseri.contigs.protein.fa \
-subject data/references/chlamy.5.3.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/reinhardtii-tblastn/C_schloesseri_hits.tsv

time tblastn \
-query data/fastas/E_debaryana.contigs.protein.fa \
-subject data/references/chlamy.5.3.fa \
-outfmt 6 \
-evalue 0.01 \
-max_target_seqs 1 \
-out data/reinhardtii-tblastn/E_debaryana_hits.tsv
```

## 5/1/2020

took between 40m-2h, with debaryana taking the longest

at a glance, seems there are hits all over the genome -
many hits seem concentrated in chrs 6 (outside the MT locus
coords) and chromosomes 13-14

need to download these files and have a look in an Rmd,
though I think it might be good to filter out contigs
where only one gene had a hit in the original reciprocal best
blast/tblastn and focus on the likely candidates for MT loci
(eg C0033 in incerta)

## 16/1/2020

alright, back on this

bringing the three files above offline to analyze in an Rmd

in the meantime - attempting to install jcvi for synteny analysis:

```bash
pip install --user jcvi
```

note - this is a python2 thing and dropped a ton of warnings

using this -

```bash
python -m jcvi.formats.fasta -h
```

apparently dropping just `python2.7 -m jcvi' raises an error - need to call
on a module each time

going to switch over to a jupyter notebook (`synteny_plots.ipynb`) in `analysis` for this

## 22/1/2020

missing some C domain genes in the original blasting! 

after looking at NameTranslation, it seems it's really just MAT3 - it just was missing
from NameTranslation, despite being in mtPlus_CDS.fasta

the PACid is 26893469 - will just pull out that sequence from the fasta
and append it to the files in `fastas-nuc` and `fastas`

```python
>>> mat3_id = '26893469'
  2 for record in SeqIO.parse('fastas-nuc/mtPlus_CDS.fasta', 'fasta'):
  3     if str(record.id) == mat3_id:
  4         mat3 = record
  5         break
>>> mat3
SeqRecord(seq=Seq('ATGTCAACCACGCACCCGCCAGAGCGCGGGCTTGTAGCACTAATTAAGGGGCTG...TGA', SingleLetterAlphabet()), id='26893469', name='26893469', description='26893469', dbxrefs=[])
>>> mat3_protein = mat3.seq.translate()
>>> with open('fastas/mtMinus_proteins.fasta', 'a') as f:
  2     f.write('>MAT3\n')
  3     f.write(str(mat3_protein) + '\n')
```

and now to run the tblastn and blastp operations again - see 28/12/2019 for code

WAIT - it seems none of the C domain genes in NameTranslation were put into mtMinus_proteins!
the reason for this is that my original name traslation worked off of the ADF IDs for mtMinus,
and no C domain gene has one of those

since the C domain genes are virtually identical between the two, we'll just pull
the CDS from mtPlus_CDS.fasta using the 268 series of IDs, and then add those sequences
to mt

finally, for the jcvi stuff, need to provide an MT only gff for reinhardtii - use mtPlus to 
keep the C domain coordinates consistent, however

prepping `mtMinus_proteins.fasta`:

```python
>>> names = {}
  2 with open('data/NameTranslation_corrected.txt', 'r') as f:
  3     for line in f:
  4         if line.startswith('C'):
  5             sp = line.split('\t')
  6             names[sp[2]] = sp[1]
>>> names
{'': '522917', '26893775': '522918', '26894611': '294742', '26893186': 'FUM1', '26894697': '522914', '26893900': 'Cre06.g255250', '26893692': 'THI10', '26893448': 'MT0829', '26893018': '522919', '26893556': '522922', '26894174': 'FBX9', '26894024': '344092', '26893348': 'MT0828', '26893097': '196063', '26893146': '522915', '26894287': 'SAD1', '26893271': '294752', '26893675': 'CGLD28', '26893528': '196073', '26893469': 'MAT3', '26894071': 'Cre06.g255150'}
>>> names.pop('')
'522917'
>>> from Bio import SeqIO
>>> from tqdm import tqdm
  2 with open('data/fastas/mtMinus_proteins.fasta', 'a') as f:
  3     for record in tqdm(SeqIO.parse(fname, 'fasta')):
  4         if str(record.id) in names.keys():
  5             f.write('>' + names[str(record.id)] + '\n')
  6             f.write(str(record.seq.translate() + '\n'))
114it [00:00, 5772.67it/s]
```

and now the blast commands can be run again.

## 24/1/2020

today - adding the C domain to mtMinus

according to Fig 1 of de Hoff, MT0828 is the first gene in the C domain,
and is in both mtPlus_only.gff (listed as gene 26893348, line 1607) and
mtMinus.gff (line 343) although there seems to be a disparity in their lengths

the mtMinus GFF only lists CDS regions - let's still line up the end coordinates
of the two gene records (345366 for MT-, 826768 for MT+) and then add the remaining
records

for jcvi with MT+, I used mRNA records with the ness_ID field in INFO - 
I ought to add gene records to the new mtMinus file, and the correct
gene names (use NameTranslation for this) in the gene field in INFO 
to be consistent with the remainder of MT-

```bash
cp -v data/references/mtMinus.gff data/references/mtMinus_C.gff
```

```python
>>> name_translation
{'': '522917', '26893775': '522918', '26893675': 'CGLD28', '26894071': 'Cre06.g255150', 
'26893271': '294752', '26893469': 'MAT3', '26894611': '294742', 
'26893900': 'Cre06.g255250', '26893448': 'MT0829', '26894174': 'FBX9', 
'26893186': 'FUM1', '26893556': '522922', '26893348': 'MT0828', 
'26893528': '196073', '26893018': '522919', '26894697': '522914', 
'26894287': 'SAD1', '26893097': '196063', '26894024': '344092', 
'26893692': 'THI10', '26893146': '522915'}
>>> name_translation.pop('') # remove one gene name w/o a ness ID
'522917'
>>> from tqdm import tqdm
  2 from copy import deepcopy
  3 import re
  4 dupes = ['26894286', '26893270']
  5 name_translation['26893507'] = '112569' # missing gene from earlier
  6 with open('data/references/mtPlus_only.gff', 'r') as f:
  7     with open('data/references/mtMinus_C.gff', 'a') as f_out:
  8         relevant_lines = []
  9         for line in tqdm(f):
 10             sp = line.split('\t')
 11             start, end = int(sp[3]), int(sp[4])
 12             if start > 826768:
 13                 if sp[2] == 'mRNA':
 14                     gene_name = sp[8].split(';')[0].lstrip('ID=PAC:')
 15                     if gene_name in dupes:
 16                         continue # some genes recorded twice for some reason
 17                     actual_name = name_translation[gene_name]
 18                     sp_to_write = deepcopy(sp)
 19                     sp_to_write[0] = 'mtMinus'
 20                     sp_to_write[1] = 'feature'
 21                     sp_to_write[2] = 'gene'
 22                     sp_to_write[3:5] = start - 481402, end - 481402
 23                     sp_to_write[8] = 'gene={};ID={}'.format(actual_name, gene_name)
 24                     out = [str(i) for i in sp_to_write]
 25                     out = '\t'.join(out)
 26                     f_out.write(out + '\n')
2168it [00:00, 400583.75it/s]
```

looks good, although this brought to my attention that we're missing _another_
gene for the mtMinus proteins/fastas - 112569, the second last gene in the C domain

```python
>>> from Bio import SeqIO
>>> with open('data/fastas/mtMinus_proteins.fasta', 'a') as f:
  2     for record in SeqIO.parse('data/fastas-nuc/mtPlus_CDS.fasta', 'fasta'):
  3         if str(record.id) == '26893507': # 112569
  4             f.write('>' + '112569' + '\n')
  5             f.write(str(record.seq.translate() + '\n'))
```

and now to blast *again* for what is hopefully the final time - once again,
the commands are up in the 28/12/2019 log

also updated NameTranslation_corrected with this missing gene

## 26/1/2020

also need to update mtMinus named fasta file for synteny plots

```python
>>> fname = 'data/NameTranslation_corrected.txt'
>>> out = 'data/fastas-nuc/mtMinus_CDS_named.fasta'
>>> plus = 'data/fastas-nuc/mtPlus_CDS.fasta'
>>> from Bio import SeqIO
>>> names = {}
  2 with open(fname, 'r') as f:
  3     for line in f:
  4         if line.startswith('C'):
  5             sp = line.split('\t')
  6             names[sp[2]] = sp[1]
>>> names.pop('')
'522917'
>>> from tqdm import tqdm
  2 with open(out, 'a') as f:
  3     for record in tqdm(SeqIO.parse(plus, 'fasta')):
  4         if str(record.id) in names.keys():
  5             f.write('>' + names[str(record.id)] + '\n')
  6             f.write(str(record.seq) + '\n')
114it [00:00, 1659.46it/s]
>>> new_names = [s.id for s in SeqIO.parse(out, 'fasta')]
>>> len(new_names)
62
>>> len(list(set(new_names)))
61
>>> for i, element in enumerate(new_names):
  2     if new_names.count(element) > 1:
  3         print(element)
MT0828
MT0828
```

will manually remove MT0828 from the file - they seem to differ slightly from
the eye test, so I'm going to keep the MT- version from Rob's files

## 27/1/2020

some manual gene changes to `mtMinus_C_noM.gff` and/or the corresponding fasta:

- HRGP1 is misspelled as HGRP1 in the mtMinus gff - correcting this manually
- changing SPP1C to SPP3 in the gff (matching the fasta + de Hoff)
- ADF43182 is actually 155027 - need to change that in the fasta (gff has the right name)
    - doing the same for `fastas/mtMinus_proteins.fasta`
- gff lists MT0618 as MTP0618 - de Hoff calls it MT0618
    - same for MT0796/MTP0796

fixed all the above manually with vim - maybe not the most automatically reproducible, but
it got the job done

rerunning the blast commands one final time after all these minor name fixes

## 13/2/2020

today: creating extended synteny diagrams

I've created a subset of chr6 for MT+ in `data/chr6_subset.gff`, but there are three
tasks remaining:

1. extract CDSs for these for the alignment
2. create a spliced MT- gff with updated coordinates
3. extract CDSs for the spliced MT- gff as well

the overall CDS file can be found in 
`/scratch/research/references/chlamy/5.3.../annotation/Creinhardtii_v5.3_223_cds.fa`

going to create a symlink in `data/fastas-nuc` and work from there

something to note - the MT- gff only contains 'ADF' format IDs - these line up
with the mtMinus fastas we have, but not the autosomes, so when creating a fasta
for jcvi to use in creating a synteny plot we have to have ADF labels for
MT- genes and standard IDs (g[0-9] etc) for the CDSs

extracting CDSs:

```python
from tqdm import tqdm
from Bio import SeqIO
seq_file = 'data/fastas-nuc/Creinhardtii_v5.3_223_cds.fa'
seqs = {}
for record in tqdm(SeqIO.parse(seq_file, 'fasta')):
    pacid = str(record.id).split('|')[1].lstrip('PACid:')
    seqs[pacid] = str(record.seq)
with open('data/fastas-nuc/chr6_subset_plus.fasta', 'w') as f_out:
    with open('data/references/chr6_subset.gff', 'r') as f:
        for line in tqdm(f):
            sp = line.rstrip().split('\t')
            chrom, type, info = sp[0], sp[2], sp[8]
            gene_id = info.split(';')[0].lstrip('ID=PAC:')
            if gene_id in seqs:
                f_out.write('>' + gene_id + '\n')
                f_out.write(seqs[gene_id] + '\n')
```
        
plot's been made in jupyter and looks great! seems regions outside MT are quite syntenic
for the most part. now to redo with mtMinus for completeness (ie steps 2 and 3 above)

## 17/2/2020

today: continuing extended synteny diagrams, but will need to 'splice in' mtMinus

comparing the plus and minus gffs (`mtPlus_genes.gff` and `mtMinus_C_noM.gff`, the final
gene of the C domain (MAT3) should have the same length in both, and it does
(chr6:937146-943474 and mtMinus:455744-462072, both 6328 bp). 

I used the ness_ID for the plus file - it's essentially the PACid minus 'PAC:'. we'll
continue to use the ness_ID, but the `.cds` file being made for jcvi should have the
mtMinus genes in 'ADF' format like in the gff

back to the offset - 943474 - 462072 = 481402, which is the value that all
post-mtMinus gene positions have to be subtracted by in the spliced gff - but also
have to account for the NIC7 offset at the start (336959 in MT+, 39983 in MT-) -
add 296976 to each MT- coord

also going to convert the type field from 'gene' to 'mRNA' to keep things consistent
(since that's what I worked off in the plus)

```python
from tqdm import tqdm
chr6_fname = 'data/references/chr6.gff'
mtminus_fname = 'data/references/mtMinus_C_noM.gff'
outname = 'data/references/chr6_mtMinus_full.gff' # will make subset later
with open(outname, 'w') as f_out:
    with open(chr6_fname, 'r') as f_chr:
        with open(mtminus_fname, 'r') as f_mt:
            counter = 0
            for line in tqdm(f_chr):
                if not 'mRNA' in line:
                    continue
                sp = line.rstrip().split('\t')
                start, end = int(sp[3]), int(sp[4])
                if start < 336959: # NIC7 start
                    f_out.write(line)
                elif start >= 336959 and end <= 943474:
                    for mt_line in f_mt:
                        sp_mt = mt_line.rstrip().split('\t')
                        if len(sp_mt[-1].split(';')) == 2:
                            sp_mt[-1] = sp_mt[-1].replace('ID=', 'ness_ID=') # C domain genes
                        m_start, m_end = int(sp_mt[3]) + 296976, int(sp_mt[4]) + 296976
                        m_start, m_end = str(m_start), str(m_end)
                        sp_out = '\t'.join(['chromosome_6', 'phytozome8_0', 'mRNA']) + '\t'
                        sp_out += '\t'.join([m_start, m_end] + sp_mt[5:])
                        sp_out += '\n'
                        f_out.write(sp_out)
                        counter += 1
                elif start > 943474:
                    sp = line.rstrip().split('\t')
                    e_start, e_end = int(sp[3]) - 481402, int(sp[4]) - 481402
                    e_start, e_end = str(e_start), str(e_end)
                    sp_out = '\t'.join(sp[:3] + [e_start, e_end] + sp[5:]) + '\n'
                    f_out.write(sp_out)
```

it works! though that took *way* longer to debug than I'd care to admit. 

creating a fasta that has all the included genes should just involve slamming `chr6_subset_plus.fasta`
together with `mtMinus_CDS.fasta` (since that has all the ADF genes with that naming):

```bash
cat chr6_subset_plus.fasta mtMinus_CDS.fasta > chr6_subset_minus.fasta
```

creating the subsetted gff (should end on the same gene as the plus - PACid 26893493) - 
corresponds to 723162 in the GFF - let's do this with awk

```bash
awk '{ if($5 <= 723162) { print }}' chr6_mtMinus_full.gff > chr6_mtMinus_subset.gff
```

back to the notebook! 

## 21/2/2020

two issues with the genes in `chr6_subset_minus.fasta` -

(well, more like 1.5)

- some of the ADF gene names have a .1 appended to the end in the bed file
- the mtPlus gene ids are going to raise 'not in `chr6_subset_minus.bed` errors
  since the MT- locus is spliced in
    - this is a much more 'ignorable' issue

fixing the ADF issue:

```python
from Bio import SeqIO
from tqdm import tqdm
import re
fname = 'data/fastas-nuc/chr6_subset_minus.fasta'
outname = 'data/fastas-nuc/chr6_subset_minus_corrected.fasta'

with open(outname, 'w') as f:
    outer_count, inner_count = 0, 0
    for record in tqdm(SeqIO.parse(fname, 'fasta')):
        seqname = str(record.id)
        pattern = '\.[0-9]$'
        if seqname.startswith('ADF'):
            outer_count += 1
            if re.search(pattern, seqname):
                seqname = seqname[:-2]
                inner_count += 1
        f.write('>' + seqname + '\n')
        f.write(str(record.seq) + '\n')
print(outer_count, inner_count) # both 41
```

# 30/4/2020

getting back on this - let's try chr6:1-4m, instead of just the first 1.2 Mb

- for `data/references/chr6.gff`, the last gene boundary is at 4004861 (line 16274 - gene 26893365)
- for `data/references/chr6_mtMinus_full.gff`, the last gene boundary is at 3993587 - should be fine

creating fastas (modified from 13/2 and 17/2 code):

```python
from tqdm import tqdm
from Bio import SeqIO
seq_file = 'data/fastas-nuc/Creinhardtii_v5.3_223_cds.fa'
seqs = {}
for record in tqdm(SeqIO.parse(seq_file, 'fasta'):
    pacid = str(record.id).split('|')[1].lstrip('PACid:')
    seqs[pacid] = str(record.seq)
with open('data/fastas-nuc/chr6_4m_plus.fasta', 'w') as f_out:
    with open('data/references/chr6.gff', 'r') as f:
        for line in tqdm(f):
            sp = line.rstrip().split('\t')
            chrom, type, end, info = sp[0], sp[2], int(sp[4]), sp[8]
            print(end)
            if end > 4004861:
                break
            gene_id = info.split(';')[0].lstrip('ID=PAC:')
            if gene_id in seqs:
                f_out.write('>' + gene_id + '\n')
                f_out.write(seqs[gene_id] + '\n') # 822 genes written to file
```


and now to create gffs that only go up to 4 Mb:

```bash
# plus
awk '{ if(($5 <= 4004861) && ($3 == "mRNA")) { print }}' data/references/chr6.gff > data/references/chr6_mtPlus_4m.gff

# minus
awk '{ if($5 <= 4000000) { print }}' data/references/chr6_mtMinus_full.gff > data/references/chr6_mtMinus_4m.gff
```

## 1/5/2020

today:

- fix dotplot axis labels
- create a minus fasta that goes up to 4m or so
- check jcvi LAST best hits against your own

from yesterday's plus_4m + incerta analysis, saw that 26893313 was the last hit at the 
3' end (matching g2240.t1) and that this region was pretty much syntenic with incerta

eyeballing `chr6_mtMinus_4m.gff`, it seems the chr6 genes here are also syntenic with the mtPlus file
(as expected, given that this is about 1 Mb in and part of the autosome) which means it's likely safe to assume
this is as far as we'll need to go for the final MT- plot (ie don't need to go all the way out to 4m)

should also do this with the schloesseri contig and see how far that goes

update: schloesseri seems to go as far as 26894644, about 20 genes upstream of 3313 - makes since given
the schloesseri contig is 0.92 Mb and the incerta contig is 1.23 Mb

I can't seem to fix these dotplot axis labels though - that said, I could
probably use these `last` and `last.filtered` files alongside the concatenated
`.bed` files to create my own plots in R (since the latter has coordinate info)

though the `last` files don't have genomic coordinates, which will have to be fetched from
the `bed` files

## 11/5/2020

today: check jcvi LAST best hits against reciprocal blastp hits

the hits files are labelled `.last.filtered` while the coords can
be found in the `.bed` files generated by jcvi

the best way to do this is to make use of `synteny_reciprocal.csv`, generated
in the R notebooks, but translate the common gene names to the
same gene IDs jcvi has (e.g. the `ID` key in the gff, instead of `gene`)

need to first split `data/synteny_reciprocal.csv` into incerta and
schloesseri files:

```bash
(head -n 1 synteny_reciprocal.csv; grep 'incerta' synteny_reciprocal.csv') > synteny_reciprocal_incerta.csv
(head -n 1 synteny_reciprocal.csv; grep 'schloesseri' synteny_reciprocal.csv') > synteny_reciprocal_schloesseri.csv
```

quick python script to append these IDs to the synteny files:

```python
from tqdm import tqdm
import csv
import re

fname = 'data/synteny_reciprocal_{}.csv'
outname = 'data/synteny_reciprocal_{}_id.csv'
gff = 'data/references/mtMinus_C_noM.gff'

# create gff lookup
gene_lookup = {}
with open(gff, 'r') as f:
    for line in f:
        info_col = line.rstrip('\n').split('\t')[-1]
        ids = info_col.split(';')
        common_name, id_name = None, None
        for gene_id in ids:
            if gene_id.startswith('gene'):
                common_name = gene_id.split('=')[1]
                if common_name.endswith('m'):
                    common_name = common_name[:-1]
                else:
                    pass
            if gene_id.startswith('ID='):
                if 'Cre' in gene_id:
                    id_name = gene.id.split('=')[1] # for Cre06 genes
                else:
                    id_name = re.search('^ID=(\w+)', gene_id).group(1)
        gene_lookup[common_name] = id_name

for species in ['incerta', 'schloesseri']:
    with open(fname.format(species), 'r', newline='') as f_in:
        reader = csv.DictReader(f_in)
        fieldnames = reader.fieldnames
        with open(outname.format(species), 'w') as f_out:
            fieldnames.append('gene_id')
            writer = csv.DictWriter(f_out, fieldnames=fieldnames)
            writer.writeheader()
            for line in tqdm(reader):
                gene_id = gene_lookup[line['gene_name']]
                line_out = line
                line_out['gene_id'] = gene_id
                writer.writerow(line_out)

```

edit: not so quick - that took a while to debug...

now to compare the two files for incerta and schloesseri. 

```R
library(tidyverse)

inc_id <- read_csv('data/synteny_reciprocal_incerta.csv', col_types = cols())
last_id <- read_tsv('jcvi/mtMinus.C_incerta.last.filtered', 
    col_types = cols(), col_names = FALSE) %>%
    select(gene_name = X1, ortholog_name = X2) %>%
    mutate(ortholog_name = str_replace(ortholog_name, '\\.t\\d', ''))

inc_id %>% select(gene_name, ortholog_name, contig) %>% anti_join(last_id) %>% as.data.frame()
```

it seems the only gene that doesn't line up is Cre06.g255150, which didn't have
a reciprocal best blast hit in incerta but showed up in the LAST alignment - interesting
too that it had reciprocal best blast hits in schloesseri and debaryana but not incerta - will
have to check the original blastp files to see if it even shows up at all. otherwise, all the genes
line up nicely! 

## 12/5/2020

turns out I didn't even need to do the name translation above - the `mtMinus.C_incerta` files
have the original gene names. whoops! 

also, HDH1 and SPP3 don't have alignments in the LAST file, but do in blastp -
otherwise, all the alignments match up cleanly. both of these line up with
contigs that aren't the 'dominant' ones - e.g. C0033 for incerta - which means
they were likely translocated elsewhere, especially considering that there
aren't any 'secondary' matches in the LAST data

Cre06.255150 is just not in the blastp alignment at all, but the match (g2152)
is on C0033 in incerta

repeating the check above with schloesseri:

```R
library(tidyverse)

sch_id <- read_csv('data/synteny_reciprocal_schloesseri.csv', col_types = cols())
last_id <- read_tsv('jcvi/mtMinus.C_schloesseri.last.filtered',
    col_types = cols(), col_names = FALSE) %>%
    select(gene_name = X1, ortholog_name = X2) %>%
    mutate(ortholog_name = str_replace(ortholog_name, '\\.t\\d', ''))

sch_id %>% select(gene_name, ortholog_name, contig) %>% anti_join(last_id) %>% as.data.frame()
```

here, SPP3 and LPS1 have matches on C0046 and C0051 respectively in the original blastp,
but do have matches in the LAST alignment as well. I'm going to trust the former, since the
latter was just done with C0045 provided, and so it's likely that these hits are 'secondary
alignments' that are only showing up since the other contigs were not present in this alignment

these translocations should be listed in the supplementary presence/absence table, so I think
they can be safely excluded from the dotplots

now to generate files to create dotplots from - this will essentially involve creating
new versions of the `last.filtered` files the the query/subject coordinates replaced with
genomic coordinates using the mtMinus and other species' GFFs

## 14/5/2020

today - generating files to create dotplots from

I think it's probably best to actually get the `last.filtered` and `bed` files offline to work with in RStudio
from the get go

```bash
mkdir data/jcvi-files
cp -v jcvi/mtMinus.C_incerta.bed data/jcvi-files/
cp -v jcvi/mtMinus.C_incerta.last.filtered data/jcvi-files/
cp -v jcvi/mtMinus.C_schloesseri.bed data/jcvi-files/
cp -v jcvi/mtMinus.C_schloesseri.last.filtered data/jcvi-files/
```

## 16/5/2020

need to repeat the above with more autosomal sequence

the tricky thing here is that the chr6 gene IDs will have changed from the
actual gene names to PACids - might need to create a name lookup table for this
just to make sure everything lines up correctly with the blastp results

```bash
cp -v jcvi/chr6_subset_minus.C_incerta.bed data/jcvi-files/
cp -v jcvi/chr6_subset_minus.C_incerta.last.filtered data/jcvi-files/
cp -v jcvi/chr6_subset_minus.C_schloesseri.bed data/jcvi-files/
cp -v jcvi/chr6_subset_minus.C_schloesseri.last.filtered data/jcvi-files/
```

update - this isn't enough! the last match in `chr6_mtPlus_4m` was around chr6:1.46m (gene 26894263) -
need to make extended chr6 subset bed and last files 

using `data/references/chr6_mtMinus_4m.gff`

need to resolve issue of no pacids in mtMinus genes - need to use the ness IDs for these

```python
from tqdm import tqdm
from Bio import SeqIO
seq_file = 'data/fastas-nuc/Creinhardtii_v5.3_223_cds.fa'
seqs = {}
for record in tqdm(SeqIO.parse(seq_file, 'fasta')):
    pacid = str(record.id).split('|')[1].lstrip('PACid:')
    seqs[pacid] = str(record.seq)
with open('data/fastas-nuc/chr6_4m_minus_temp.fasta', 'w') as f_out:
    with open('data/references/chr6_mtMinus_full.gff', 'r') as f:
        for line in tqdm(f):
            if 'ID=PAC:' not in line:
                continue # will append mtMinus_CDS.fasta to get these
            sp = line.rstrip().split('\t')
            chrom, type, end, info = sp[0], sp[2], int(sp[4]), sp[8]
            if end > 3993587:
                print(end)
                break
            gene_id = info.split(';')[0].lstrip('ID=PAC:')
            if gene_id in seqs:
                f_out.write('>' + gene_id + '\n')
                f_out.write(seqs[gene_id] + '\n')
```

and then:

```bash
# in data/fastas-nuc
cat chr6_4m_minus_temp.fasta mtMinus_CDS.fasta > chr6_4m_minus.fasta
rm chr6_4m_minus_temp.fasta
grep -c '>' chr6_4m_minus.fasta # 834 genes
```

now to run the jcvi commands:

```bash
cd jcvi/

python2.7 -m jcvi.formats.gff bed --type=mRNA \
--key=ness_ID ../data/references/chr6_mtMinus_4m.gff \
-o chr6_mtMinus_4m_temp.bed

python2.7 -m jcvi.formats.fasta format \
../data/fastas-nuc/chr6_4m_minus.fasta \
chr6_mtMinus_4m.cds
```

need to remove the .1s from the ADF genes listed in the bed file:

```python
import re
from copy import deepcopy
from tqdm import tqdm
fname = 'jcvi/chr6_mtMinus_4m_temp.bed'
outname = 'jcvi/chr6_mtMinus_4m.bed'
from copy import deepcopy
from tqdm import tqdm
import re
with open(fname, 'r') as f_in:
    fieldnames = ['chrom', 'start', 'end', 'gene', 'score', 'strand']
    reader = csv.DictReader(f_in, fieldnames=fieldnames, delimiter='\t')
    with open(outname, 'w') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        for line in tqdm(reader):
            line_out = deepcopy(line)
            if 'ADF' in line['gene']:
                line_out['gene'] = re.search('(ADF[0-9]+)', line['gene']).group(1)
                writer.writerow(line_out)
            else:
                writer.writerow(line_out)
```

and now we can do the alignments for incerta and schloesseri:

```bash
time python2.7 -m jcvi.compara.catalog ortholog \
chr6_mtMinus_4m C_incerta --cscore=.99

time python2.7 -m jcvi.compara.catalog ortholog \
chr6_mtMinus_4m C_schloesseri --cscore=.99
```

porting the files to `data/jcvi-files`:

```bash
cat chr6_mtMinus_4m.bed C_incerta.bed > chr6_mtMinus_4m.C_incerta.bed
cat chr6_mtMinus_4m.bed C_schloesseri.bed > chr6_mtMinus_4m.C_schloesseri.bed

cp -v chr6_mtMinus_4m.C_incerta.bed ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_incerta.last.filtered ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_schloesseri.bed ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_schloesseri.last.filtered ../data/jcvi-files
```

## 25/5/2020

I goofed! I goofed so hard making `chr6_mtMinus_full.gff` (see 17/2 log)

the offset is incorrect - instead of doing 943k - 462k, I have to start the gene
after MAT3 at the same distance (?) downstream _relative_ to MAT3 that it is
in the MT+ file, which the previous method did not do - leading to the gene
immediately after it being at like 462k instead of being greater than
the MAT3 end coord of ~759k

the *actual* offset post-MAT3 is 943474 - (759048 + 3) = 184423, since the next gene
starts just 3 bp downstream at 943477 in MT+

this code should handle the offset correctly:

```python
from tqdm import tqdm
chr6_fname = 'data/references/chr6.gff'
mtminus_fname = 'data/references/mtMinus_C_noM.gff'
outname = 'data/references/chr6_mtMinus_full.gff' # will make subset later
with open(outname, 'w') as f_out:
    with open(chr6_fname, 'r') as f_chr:
        with open(mtminus_fname, 'r') as f_mt:
            counter = 0
            for line in tqdm(f_chr):
                if not 'mRNA' in line:
                    continue
                sp = line.rstrip().split('\t')
                start, end = int(sp[3]), int(sp[4])
                if start < 336959: # NIC7 start
                    f_out.write(line)
                elif start >= 336959 and end <= 943474:
                    for mt_line in f_mt:
                        sp_mt = mt_line.rstrip().split('\t')
                        if len(sp_mt[-1].split(';')) == 2:
                            sp_mt[-1] = sp_mt[-1].replace('ID=', 'ness_ID=') # C domain genes
                        m_start, m_end = int(sp_mt[3]) + 296976, int(sp_mt[4]) + 296976
                        m_start, m_end = str(m_start), str(m_end)
                        sp_out = '\t'.join(['chromosome_6', 'phytozome8_0', 'mRNA']) + '\t'
                        sp_out += '\t'.join([m_start, m_end] + sp_mt[5:])
                        sp_out += '\n'
                        f_out.write(sp_out)
                        counter += 1
                elif start > 943474:
                    sp = line.rstrip().split('\t')
                    e_start, e_end = int(sp[3]) - 184423, int(sp[4]) - 184423
                    e_start, e_end = str(e_start), str(e_end)
                    sp_out = '\t'.join(sp[:3] + [e_start, e_end] + sp[5:]) + '\n'
                    f_out.write(sp_out)

```

remaking the gff:

```bash
awk '{ if($5 <= 4000000) { print }}' data/references/chr6_mtMinus_full.gff > data/references/chr6_mtMinus_4m.gff
```

rerunning 16/5 commands:

```python
from tqdm import tqdm
from Bio import SeqIO
seq_file = 'data/fastas-nuc/Creinhardtii_v5.3_223_cds.fa'
seqs = {}
for record in tqdm(SeqIO.parse(seq_file, 'fasta')):
    pacid = str(record.id).split('|')[1].lstrip('PACid:')
    seqs[pacid] = str(record.seq)
with open('data/fastas-nuc/chr6_4m_minus_temp.fasta', 'w') as f_out:
    with open('data/references/chr6_mtMinus_full.gff', 'r') as f:
        for line in tqdm(f):
            if 'ID=PAC:' not in line:
                continue # will append mtMinus_CDS.fasta to get these
            sp = line.rstrip().split('\t')
            chrom, type, end, info = sp[0], sp[2], int(sp[4]), sp[8]
            if end > 3993587:
                print(end)
                break
            gene_id = info.split(';')[0].lstrip('ID=PAC:')
            if gene_id in seqs:
                f_out.write('>' + gene_id + '\n')
                f_out.write(seqs[gene_id] + '\n')

```

and then:

```bash
cd data/fastas-nuc
cat chr6_4m_minus_temp.fasta mtMinus_CDS.fasta > chr6_4m_minus.fasta
rm chr6_4m_minus_temp.fasta
grep -c '>' chr6_4m_minus.fasta # 781 genes! 

cd ../../jcvi
python2.7 -m jcvi.formats.gff bed --type=mRNA \
--key=ness_ID ../data/references/chr6_mtMinus_4m.gff \
-o chr6_mtMinus_4m_temp.bed

python2.7 -m jcvi.formats.fasta format \
../data/fastas-nuc/chr6_4m_minus.fasta \
chr6_mtMinus_4m.cds
```

cleaning the gene names:

```python
import re
from copy import deepcopy
from tqdm import tqdm
import csv
fname = 'jcvi/chr6_mtMinus_4m_temp.bed'
outname = 'jcvi/chr6_mtMinus_4m.bed'
with open(fname, 'r') as f_in:
    fieldnames = ['chrom', 'start', 'end', 'gene', 'score', 'strand']
    reader = csv.DictReader(f_in, fieldnames=fieldnames, delimiter='\t')
    with open(outname, 'w') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        for line in tqdm(reader):
            line_out = deepcopy(line)
            if 'ADF' in line['gene']:
                line_out['gene'] = re.search('(ADF[0-9]+)', line['gene']).group(1)
                writer.writerow(line_out)
            else:
                writer.writerow(line_out)
```

redoing alignments:


```bash
time python2.7 -m jcvi.compara.catalog ortholog \
chr6_mtMinus_4m C_incerta --cscore=.99

time python2.7 -m jcvi.compara.catalog ortholog \
chr6_mtMinus_4m C_schloesseri --cscore=.99
```

replacing the files in `data/jcvi-files`:

```bash
rm data/jcvi-files/chr6_mtMinus_4m.C_*

cd jcvi/
cat chr6_mtMinus_4m.bed C_incerta.bed > chr6_mtMinus_4m.C_incerta.bed
cat chr6_mtMinus_4m.bed C_schloesseri.bed > chr6_mtMinus_4m.C_schloesseri.bed

cp -v chr6_mtMinus_4m.C_incerta.bed ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_incerta.last.filtered ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_schloesseri.bed ../data/jcvi-files
cp -v chr6_mtMinus_4m.C_schloesseri.last.filtered ../data/jcvi-files
```










