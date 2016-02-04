#!/bin/bash

## COMMANDS FOR DOWNLOADING ADDITIONAL REQUIRED REFERENCE FILES FOR CONVERTING CNV GENOTYPES IN GENOMESTRIP VCF TO PLINK FORMAT.

## For Drosophila_melanogaster build dm6. Requires liftover of dm3 files.

## Primarily for Linux use but mostly transferable to mac/unix.

## see http://parsecnv.sourceforge.net/
mkdir GeneRef_DM6
cd GeneRef_DM6

wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/refLink.txt.gz
###wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/gc5BaseBw.txt.gz
wget http://hgdownload.cse.ucsc.edu/gbdb/dm6/bbi/gc5BaseBw/gc5Base.bw
###Use hg19 Pathways files by querying TOUPPER of all genes
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/cytoBandIdeo.txt.gz

wget http://humanparalogy.gs.washington.edu/dm3/data/dm3genomicSuperDup.tab ###Needs Liftover as Instructed Below ### No Header See Header in wget http://paralogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
wget http://paralogy.gs.washington.edu/build37/data/GRCh37GenomicSuperDup.tab
head -1 GRCh37GenomicSuperDup.tab > header_SuperDub.tab
cat header_SuperDub.tab dm3genomicSuperDup.tab > dm3genomicSuperDup_headed.tab
mv dm3genomicSuperDup_headed.tab dm3genomicSuperDup.tab


wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Drosophila_melanogaster/by_assembly/Release_6_plus_ISO1_MT/gvf/Release_6_plus_ISO1_MT.remap.all.germline.ucsc.gvf.gz
gunzip refGene.txt.gz; mv refGene.txt dm6_refGene.txt
gunzip refLink.txt.gz; mv refLink.txt dm6_refLink.txt
###gunzip gc5BaseBw.txt.gz; mv gc5BaseBw.txt dm6_gc5BaseBw.txt ## did this.

cp ../PolishRefGeneForScanRegion_KeepOnlyChrLines.pl .
perl PolishRefGeneForScanRegion_KeepOnlyChrLines.pl dm6_refGene.txt ### to delete non “chr” lines
cp ../GeneRef/MakeExonScanRegionDefFile* .
perl MakeExonScanRegionDefFile_refGene.pl dm6_refGene.txt_OnlyChr.txt
perl MakeExonScanRegionDefFile_wGene_refGene.pl dm6_refGene.txt_OnlyChr.txt

#cp ../PerlModules/scan_region.pl . ## did this.
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod 777 bigWigToBedGraph
mv gc5Base.bw dm6_gc5Base.bw
./bigWigToBedGraph dm6_gc5Base.bw dm6_gc5Base.bedGraph ### this command doesn't work on mac, presumably as executable was made under linux.
cp ../Polishgc5BaseForScanRegion_KeepOnlyChrLines_5120Increments.pl .
perl Polishgc5BaseForScanRegion_KeepOnlyChrLines_5120Increments.pl dm6_gc5Base.bedGraph
### Illegal division by zero at Polishgc5BaseForScanRegion_KeepOnlyChrLines_5120Increments.pl line 30, <SAMPSHEET> line 358024431.# so no Y

awk -F"\t" '{print $4"\t"$1"\t+\t"$2"\t"$3"\t \t \t \t \t \t \t "}' dm6_gc5Base.bedGraph_OnlyChr_5120Intervals.txt > dm6_gc5Base_SimFormat_AllCol.sorted ### Now unique values, rand() was added before to $4 GC to allow duplicates in scan region since original 5 base pair interval values are only 20,40,60,80
###perl Polishgc5BaseForScanRegion_KeepOnlyChrLines.pl dm6_gc5Base_SimFormat_AllCol.sorted_randDec
rm dm6_gc5Base.bw; rm dm6_gc5Base.bedGraph; rm dm6_gc5Base.bedGraph_OnlyChr_5120Intervals.txt




gunzip  cytoBandIdeo.txt.gz; mv cytoBandIdeo.txt dm6_cytoBandIdeo.txt
awk -F"\t" '{print $1$4"\t"$1"\t+\t"$2"\t"$3"\t \t \t \t \t \t \t "}' dm6_cytoBandIdeo.txt | sed 's/chr//' > dm6_cytoBand_SimFormat_AllCol.sorted

awk '{print $1"\t"$2"\t"$3"\t"$6"#"$7"#"$8}' dm3genomicSuperDup.tab > dm3genomicSuperDup.tab_ForLiftover.txt




### Go to https://genome.ucsc.edu/cgi-bin/hgLiftOver 
## Original Assembly: dm3 New Assembly: dm6 
## Click Choose File 
## Click Upload File 
## Click View Conversions Copy File into GeneRef
mv hglft_genome_*.bed dm6_genomicSuperDup.tab

## Sort By Chr Pos. Note unix may have problems with \t as tab symbol depending on version.
sort -k 1,1 -k 2,2n dm6_genomicSuperDup.tab | sed 's/#/\t/g' > dm6_genomicSuperDup.tab.sorted 
# sort -k 1,1 -k 2,2n dm6_genomicSuperDup.tab | sed -e 's|\#|\t|g' > dm6_genomicSuperDup.tab.sorted 
# sort -k 1,1 -k 2,2n dm6_genomicSuperDup.tab | sed -e 's|#| |g' > dm6_genomicSuperDup.tab.sorted


awk '{print $5":"$6"\t"$1"\t"$4"\t"$2"\t"$3"\t \t \t \t \t \t \t "}' dm6_genomicSuperDup.tab.sorted > dm6_genomicSuperDups_SimFormat_AllCol.sorted

gunzip Release_6_plus_ISO1_MT.remap.all.germline.ucsc.gvf.gz
mv Release_6_plus_ISO1_MT.remap.all.germline.ucsc.gvf dm6_Release_6_plus_ISO1_MT.remap.all.germline.ucsc.gvf

## remove description rows. Remove "ID=". Replace ";Name" with tab. Replace "copy_number_" with "CN".
grep -v ^# dm6_Release_6_plus_ISO1_MT.remap.all.germline.ucsc.gvf | sed 's/ID=//' | sed 's/;Name/\t/' | sed 's/copy_number_/CN/' | sed 's/variation/var/' | awk '{print $9":"$3"\t"$1"\t+\t"$4"\t"$5"\t \t \t \t \t \t \t "}' > dm6_dgv_SimFormat_AllCol.sorted
