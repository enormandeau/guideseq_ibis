# GUIDEseq IBIS

GUIDE-seq analysis pipeline

## TODO

- Filter better (remove reads with multiple dsODN...)
- Remove PCR duplicates from SAM files
- Create simplified gene annotation file with extensive ranges
- Query exons instead of genes?


## Steps

- Prepare config file(s)
- Index genome
- Trim reads
- Extract UMI infos + 8 first nucleotides and fasta sequences
- Align to genome
- Remove PCR duplicates
