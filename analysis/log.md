
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

