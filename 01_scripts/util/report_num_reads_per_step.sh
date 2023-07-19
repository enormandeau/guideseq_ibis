#!/bin/bash
# Report number of reads per step

export FOLDER="04_data"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo -e "Sample\t""$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.gz | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo -e "$i\t"$(gunzip -c "$i"*.gz | wc -l | awk '{print $1/4}')
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="05_trimmed"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo "$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.gz | cut -d "_" -f 1 | sort -uV |
    while read i
    do
        echo $(gunzip -c "$i"*.gz | wc -l | awk '{print $1/4}')
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="06_extracted"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo "$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.gz | cut -d "." -f 1 | sort -uV |
    while read i
    do
        echo $(gunzip -c "$i"*.gz | grep -c "^>")
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="07_aligned"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo "$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.sam | cut -d "." -f 1 | sort -uV |
    while read i
    do
        echo $(wc -l "$i"*.sam | awk '{print $1}')
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="08_deduplicated"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo "$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.dedup | cut -d "." -f 1 | sort -uV |
    while read i
    do
        echo $(wc -l "$i"*.dedup | awk '{print $1}')
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

export FOLDER="09_sites"
echo "  Examining $FOLDER"
cd "$FOLDER"
echo "$FOLDER" > ../10_read_dropout/"$FOLDER"_reads.tsv
ls -1 *.dsb | cut -d "." -f 1 | sort -uV |
    while read i
    do
        echo $(grep -v Sample "$i"*.dsb | awk '{s = s+($6+$7)}END{print s}')
    done >> ../10_read_dropout/"$FOLDER"_reads.tsv
cd ..

paste 10_read_dropout/*.tsv > read_dropout.tsv
