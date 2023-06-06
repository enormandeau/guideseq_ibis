#!/bin/bash
# Report number of reads per step

export FOLDER="04_data"
cd "$FOLDER"
ls -1 *.gz | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(gunzip -c "$i"*.gz | wc -l | awk '{print $1/4}')
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="05_trimmed"
cd "$FOLDER"
ls -1 *.gz | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(gunzip -c "$i"*.gz | wc -l | awk '{print $1/4}')
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="06_extracted"
cd "$FOLDER"
ls -1 *.gz | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(gunzip -c "$i"*.gz | grep -c "^>")
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="07_aligned"
cd "$FOLDER"
ls -1 *.sam | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(wc -l "$i"*.sam | awk '{print $1}')
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="08_deduplicated"
cd "$FOLDER"
ls -1 *.dedup | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(wc -l "$i"*.dedup | awk '{print $1}')
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="09_sites"
cd "$FOLDER"
ls -1 *.dsb | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(grep -v Sample "$i"*.dsb | awk '{s = s+($6+$7)}END{print s}')
    done > ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

paste 10_read_dropout/*.tsv > read_dropout.tsv

cat 09_sites/* | head -1 > guideseq_ibis_report.tsv
for i in 09_sites/*.dsb
do
    grep -v "^Sample" "$i"
done >> guideseq_ibis_report.tsv
