
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







