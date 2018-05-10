#!/bin/bash

# Home of the programs
################################################################
CSAHDIR=/home/gaspari/csah_pipeline

# Filenames
################################################################
SPROT_FTC_OUT=uniprot_sprot._formatted_output.txt
TREMBL_FTC_OUT=uniprot_trembl._formatted_output.txt

SPROT_S4C_OUT=uniprot_sprot.s4c
TREMBL_S4C_OUT=uniprot_trembl.s4c

SPROT_FASTA_GZ=uniprot_sprot.fasta.gz
TREMBL_FASTA_GZ=uniprot_trembl.fasta.gz

# Consensus
SPROT_MASKEDFASTA=uniprot_sprot_masked.fst
SPROT_CSAHFASTA=uniprot_sprot_csah.fst
SPROT_CSAH=uniprot_sprot_csah.txt
SPROT_CSAHDAT=uniprot_sprot_csah.dat
TREMBL_MASKEDFASTA=uniprot_trembl_masked.fst
TREMBL_CSAHFASTA=uniprot_trembl_csah.fst
TREMBL_CSAH=uniprot_trembl_csah.txt
TREMBL_CSAHDAT=uniprot_trembl_csah.dat

# FT_CHARGE only
SPROT_FTC_MASKEDFASTA=uniprot_sprot_ftc_masked.fst
SPROT_FTC_CSAHFASTA=uniprot_sprot_ftc_csah.fst
SPROT_FTC_CSAH=uniprot_sprot_ftc_csah.txt
SPROT_FTC_CSAHDAT=uniprot_sprot_ftc_csah.dat
TREMBL_FTC_MASKEDFASTA=uniprot_trembl_ftc_masked.fst
TREMBL_FTC_CSAHFASTA=uniprot_trembl_ftc_csah.fst
TREMBL_FTC_CSAH=uniprot_trembl_ftc_csah.txt
TREMBL_FTC_CSAHDAT=uniprot_trembl_ftc_csah.dat



# Generate fasta file with sequences with potential CSAHs only
################################################################
grep "sp" $SPROT_FTC_OUT | awk -F"|" '{print $2}' | sort -u > sprot_ftc.ids
zcat $SPROT_FASTA_GZ | $CSAHDIR/getseqsbyidfromfst.pl -f sprot_ftc.ids > uniprot_sprot_ftc.fst

grep "tr" $TREMBL_FTC_OUT | awk -F"|" '{print $2}' | sort -u > trembl_ftc.ids
zcat $TREMBL_FASTA_GZ | $CSAHDIR/getseqsbyidfromfst.pl -f trembl_ftc.ids > uniprot_trembl_ftc.fst

# Run scan4csah.pl
################################################################
zcat $SPROT_FASTA_GZ  | $CSAHDIR/scan4csah.pl > $SPROT_S4C_OUT
zcat $TREMBL_FASTA_GZ | $CSAHDIR/scan4csah.pl > $TREMBL_S4C_OUT


# Getting consensus
################################################################
$CSAHDIR/csahdetect.pl --infasta=uniprot_sprot_ftc.fst --s4coutfile=uniprot_sprot.s4c --ftcoutfile=$SPROT_FTC_OUT --maskedfasta=$SPROT_MASKEDFASTA --csahfasta=$SPROT_CSAHFASTA > $SPROT_CSAH

$CSAHDIR/csahdetect.pl --infasta=uniprot_trembl_ftc.fst --s4coutfile=uniprot_trembl.s4c --ftcoutfile=$TREMBL_FTC_OUT --maskedfasta=$TREMBL_MASKEDFASTA --csahfasta=$TREMBL_CSAHFASTA > $TREMBL_CSAH

# Getting FT_CHARGE-only results 
################################################################
$CSAHDIR/csahdetect.pl --mode=F --infasta=uniprot_sprot_ftc.fst --ftcoutfile=$SPROT_FTC_OUT --maskedfasta=$SPROT_FTC_MASKEDFASTA --csahfasta=$SPROT_FTC_CSAHFASTA > $SPROT_FTC_CSAH

$CSAHDIR/csahdetect.pl --mode=F --infasta=uniprot_trembl_ftc.fst --ftcoutfile=$TREMBL_FTC_OUT --maskedfasta=$TREMBL_FTC_MASKEDFASTA --csahfasta=$TREMBL_FTC_CSAHFASTA > $TREMBL_FTC_CSAH

# Generating uniprot dat file with full annotations plus FT CSAH lines
################################################################
$CSAHDIR/retrieve_uniprot_dat_for_id_list.pl sprot_ftc.ids > sprot_ftc.dat
$CSAHDIR/insertCSAH2dat.pl -s -l $SPROT_CSAH < sprot_ftc.dat > $SPROT_CSAHDAT

# trembl data should be downloaded in pieces as the full list is too large
$CSAHDIR/gettrembldatinpieces.pl < trembl_ftc.ids > trembl_ftc.dat
$CSAHDIR/insertCSAH2dat.pl -s -l $TREMBL_CSAH < trembl_ftc.dat > $TREMBL_CSAHDAT

# Preparing files for download on website
#########################################
rm -rf ft_charge
mkdir ft_charge
mv $SPROT_FTC_CSAH $SPROT_FTC_CSAHFASTA $SPROT_FTC_MASKEDFASTA $TREMBL_FTC_CSAH $TREMBL_FTC_MASKEDFASTA $TREMBL_FTC_CSAHFASTA ft_charge/
zip -r CSAHDB.zip CSAHDB_README.TXT $SPROT_CSAHDAT $TREMBL_CSAHDAT $SPROT_CSAH $TREMBL_CSAH $SPROT_MASKEDFASTA $TREMBL_MASKEDFASTA $SPROT_CSAHFASTA $TREMBL_CSAHFASTA ft_charge/*


