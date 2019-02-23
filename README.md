# bifo-scripts
Bioinformatic scripts

Scripts include:
export-embl.pl

merges table files from Verdant, user annotation and gene defiition files to create an EMBL formatted annotation file for a chloroplast with the accessory gene_names.in file.

map-reads2ref.pl

download read files from NCBI's SRA, unpack them and map to a reference plastome with samtools. Merge multiple mapped bam files if needed, extract read depth at each reference position then export the log10 converted file ready for import to a drawing application.

