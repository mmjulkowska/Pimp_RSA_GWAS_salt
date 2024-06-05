
## Figuring out the input files:

I have had multiple files for SNP calling - but I dont remember which ones would be most appropriate to use - we would like to use as many SNPs as possible.

```bash
grep -v "^#" Pimp_CHR01.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#342896
grep -v "^#" Pimp_CHR02.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#239026
grep -v "^#" Pimp_CHR03.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#310298
grep -v "^#" Pimp_CHR04.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#238151
grep -v "^#" Pimp_CHR05.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#280313
grep -v "^#" Pimp_CHR06.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#195414
grep -v "^#" Pimp_CHR07.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
252140
grep -v "^#" Pimp_CHR08.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
263072
grep -v "^#" Pimp_CHR09.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#360134
grep -v "^#" Pimp_CHR10.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#355192
grep -v "^#" Pimp_CHR11.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#272019
grep -v "^#" Pimp_CHR12.fin16.vcf.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
#419212 
```

So in total we have 3.5 million SNPs in all of these files. 

```bash
for i in $(seq -w 1 12); do
grep -v "^#" Pimp_CHR$i.recode.vcf.result.vcf | grep -v "#" | wc -l
done

for i in $(seq -w 1 12); do
grep -v "^#" Pimp_CHR$i.MAC3.recode.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
done

for i in $(seq -w 1 12); do
grep -v "^#" Pimp_CHR$i.HM.vcf.mod.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
done

for i in $(seq -w 1 12); do
grep -v "^#" Pimp_CHR$i.HM.vcf | awk '$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $0}' | wc -l
done
```

OK - so the files Pimp_CHRXX.HM.vcf.mod.vcf	Pimp_CHRXX.HM.vcf	Pimp_CHRXX.fin16.vcf.recode.vcf are potentially good to be used for the SNPs in GWAS. Let's see what's their genotype coding, so we have (preferably) M- identifiers for GWAS (KAUST.ID)

```bash
grep '^#CHROM' Pimp_CHR01.recode.vcf | awk '{for(i=10; i<=NF; i++) print $i}'
```

We need to rename the accessions first:

```bash
bcftools query -l Pimp_CHR01.recode.vcf > original_samples.txt
touch new_samples.txt
echo "LA1868
LA1256
LA1572
LA1381
LA1577
LA1594
LA1593
LA1591
LA1359
LA1382
LA1380
LA2854
LA1597
LA1598
LA1599
LA1600
LA1601
LA1604
LA1258
LA1683
LA1728
LA1729
LA1645
LA2832
LA2836
LA1993
LB1258
LA2000
LA1374
LA0391
LC1258
LA0400
LA1349
LA1384
LA0480
LA1472
LA1561
LA1470
LA1259
LA1478
LA1589
LA1590
LA2628
LA2645
LA2647
LA2659
LA2718
LA2725
LA2804
LA0369
LA0412
LA0417
LA0420
LA0442
LA1575
LB1572
LA1576
LA1242
LA0859
LA0121
LA0443
LA0722
LA0753
LA1243
LA1245
LA1521
LA1520
LA1514
LA1585
LA1583
LA0114
LA1581
LB1581
LA1580
LA1578
LA1595
LA1579
LA1260
LA1263
LA1279
LA0100
LA1280
LB1280
LA1301
LA1332
LA1335
LA1342
LA1343
LB1343
LA1344
LA1345
LA1867
LA1742
LA1348
LA2537
LA2536
LA2533
LA2429
LA2426
LA2425
LA2423
LA2412
LA1828
LA1660
LA1661
LA1670
LA1676
LA1678
LA1679
LA1680
LA0373
LA0375
LA0376
LA1831
LA0381
LA1605
LA1607
LA1608
LA1610
LA1611
LA1612
LA1613
LA1838
LA1839
LA1841
LA1845
LA1687
LA1697
LA1719
LA1866
LA1720
LA1634
LA1428
LA1371
LA1416
LA1370
LA1628
LA1630
LA1633
LA1864
LA1636
LA1637
LA1571
LA1586
LA2539
LA2540
LA2544
LA2545
LA2546
LA2547
LA1852
LA2576
LA2578
LA2805
LA2839
LA2852
LA2857
LA2866
LA1587
LA1588
LA1992
LA1847
LA1870
LA1874
LA1921
LA1923
LA1924
LA1950
LA1987
LA2186
LA2188
LA2345
LA2347
LA2389
LA2398
LA2187
LA2903
LA2904
LA2914
LB2914
LA2915
LA2934
LA1246
LA2966
LA2974
LA2982
LA2983
LA3158
LA3159
LA3160
LA3161
LA3859
LA4138
LA1248
LA1685
LA2149
LA2170
LA2173
LA2176
LA2181
LA2184
LA2178
LA2180
LB1576
LA0397
LA0398
LA0411
LA0418
LA1237
LA1257
LA1261
LA1355
LA1357
LA1375
LA1383
LA1469
LA1562
LA1573
LB1593
LA1603
LA1614
LA1629
LA1631
LB1633
LB1634
LA1635
LA1651
LA1659
LA1682
LA1684
LB1685
LA1686
LA1689
LA1690
LA1810
LA1825
LA1835
LA1836
LA1842
LA1858
LA1861
LA1863
LA1925
LA1933
LA1936
LA2001
LA2096
LA2145
LA2146
LA2147
LA2183
LA2335
LA2346
LA2348
LA2401
LA2543
LA2646
LA2649
LA2653
LA2831
LA2833
LB2839
LA2840
LA2850
LA2851" > new_samples.txt

bcftools reheader -s new_samples.txt -o CHR01.vcf Pimp_CHR01.recode.vcf
bcftools reheader -s new_samples.txt -o CHR02.vcf Pimp_CHR02.recode.vcf
bcftools reheader -s new_samples.txt -o CHR03.vcf Pimp_CHR03.recode.vcf
bcftools reheader -s new_samples.txt -o CHR04.vcf Pimp_CHR04.recode.vcf
bcftools reheader -s new_samples.txt -o CHR05.vcf Pimp_CHR05.recode.vcf
bcftools reheader -s new_samples.txt -o CHR06.vcf Pimp_CHR06.recode.vcf
bcftools reheader -s new_samples.txt -o CHR07.vcf Pimp_CHR07.recode.vcf
bcftools reheader -s new_samples.txt -o CHR08.vcf Pimp_CHR08.recode.vcf
bcftools reheader -s new_samples.txt -o CHR09.vcf Pimp_CHR09.recode.vcf
bcftools reheader -s new_samples.txt -o CHR10.vcf Pimp_CHR10.recode.vcf
bcftools reheader -s new_samples.txt -o CHR11.vcf Pimp_CHR11.recode.vcf
bcftools reheader -s new_samples.txt -o CHR12.vcf Pimp_CHR12.recode.vcf
```

YAY! It worked! So now we have the files with correctly named genotypes!!! 

OK - it seems that fin16 file has more genotypes. Let's see what we have in RSA files. seems like fin16 file contains much more genotypes - and some of the duplicated accessions. So I guess we would be going to use HM.vcf.mod.vcf files.

Maybe I should put all of the other files into a different directory
```bash
mkdir old_SNPs
mv *.HM.hmp.txt old_SNPs/
mv *.HM.vcf old_SNPs/
mv *.MAC3.recode.vcf old_SNPs/
mv *.fin* old_SNPs/
mv *MAC3* old_SNPs/
mv *.result.vcf old_SNPs/
mv *.recode.vcf old_SNPs/
mv Pimp_* old_SNPs/
```

OK - so now we are working with the files called Pimp_CHRXX.HM.vcf.mod.vcf

## GWAS

Let's make sure that we have all of the libraries installed as necessary

```bash
# Install PLINK
conda install conda-forge::mamba
mamba install bioconda::plink
#prepare binary plink files
plink --vcf CHR01.vcf --recode --make-bed --out CHR01
```

Now - let's install other things:

```bash
# Working directory
pwd
/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo

wget https://pmglab.top/gec/data/archive/v0.2/gecV0.2.zip
unzip gecV0.2.zip
```

OK - first we need to re-calculate chromosome files into binary files:

```bash
# Check for one chromosome if it works well
java -jar /home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo/gec/gec.jar --effect-number --maf 0.05 --plink-binary CHR01 mv gec.sum CHR01.gec.sum
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

# Re-calculate plink binary files for all other chromosomes:
for i in $(seq -w 2 12); do
    plink --vcf CHR$i.vcf --recode --make-bed --out CHR$i
done
```

Now - let's calculate kinship and PCA on the effective number of SNPs for all chrommosomes:

```bash
for i in $(seq -w 1 12); do
    java -jar /home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo/gec/gec.jar --effect-number --maf 0.05 --plink-binary CHR$i mv gec.sum CHR$i.gec.sum
done

# Calculate the effective number of SNPs of all chromosomes
cat *gec.sum | grep -v Effective_Number | awk '{ sum += $2 } END { print sum }' 

#Let's combine all of the individual VCFs into one:
bgzip CHR01.vcf
bgzip CHR02.vcf
bgzip CHR03.vcf
bgzip CHR04.vcf
bgzip CHR05.vcf
bgzip CHR06.vcf
bgzip CHR07.vcf
bgzip CHR08.vcf
bgzip CHR09.vcf
bgzip CHR10.vcf
bgzip CHR11.vcf
bgzip CHR12.vcf

bcftools concat *.vcf.gz -o ALLCHR.vcf.gz -O z

gunzip *.vcf.gz


# LD prune
plink --vcf ALLCHR.vcf --keep <(cut -f1,2 trait.MR.GRC.txt) --indep-pairwise 50 10 0.2 --allow-extra-chr --double-id --out ALLCHR.filt
# Extract the prunned SNPs
plink --vcf ALLCHR.vcf --keep <(cut -f1,2 trait.MR.GRC.txt) --extract ALLCHR.filt.prune.in --recode --double-id --out ALLCHR.filt.prune

# Calculate PCA
plink --file ALLCHR.filt.prune --pca --out ALLCHR.filt.prune --allow-extra-chr

# Remove unnecessary files
rm *log *nosex

# Prepare the PCA file as cofactor
awk '{print $1,$2,1,$3,$4,$5,$6,$7}' ALLCHR.filt.prune.eigenvec > ALLCHR.filt.prune.PCA.txt

# Transpose the ped file
for i in $(seq -w 1 12); do
    plink --file FCHR$i --recode12 --output-missing-genotype 0 --transpose --out FCHR$i &
done

# Calculate kinship
plink --file ALLCHR.filt.prune --recode12 --output-missing-genotype 0 --transpose --allow-extra-chr --out ALLCHR.filt.prune


# EMMAX
#Install EMMAX
wget https://csg.sph.umich.edu//kang/emmax/download/emmax-intel-binary-20120210.tar.gz
gunzip emmax-intel-binary-20120210.tar.gz
tar -xvf emmax-intel-binary-20120210.tar


# Estimate kinship
./emmax-kin-intel64 -v -d 10 ALLCHR.filt.prune
```


Great - now we are ready for calculating the actual GWAS:

```bash 
# Let's try one trait for starters
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.MR.GRC.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o trait.MR.GRC.Chr$i &
done
```

OK - if we want informative output from our GWAS - we need to include Chr/pos info into 3rd collumn of VCF:

```bash
for i in $(seq -w 1 12); do
    awk 'BEGIN {FS=OFS="\t"} {if($1 ~ /^#/){print $0} else {$3 = $1 "_" $2; print $0}}' CHR$i.vcf > FCHR$i.vcf
done
```

Great - now ;et's do GWAS for all of the Growth Rate files:

```bash
for e in LRnG.STI LRno.GRS MR.GRS aLRG.STI aLRL.GRS LRno.GRC MR.GRC MRG.STI aLRL.GRC; do
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.$e.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o $e.Chr$i &
done
done
```

Add coordinates of the CHR_POS into ped file:

```bash
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR01.vcf > FCHR01.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR02.vcf > FCHR02.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR03.vcf > FCHR03.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR04.vcf > FCHR04.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR05.vcf > FCHR05.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR06.vcf > FCHR06.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR07.vcf > FCHR07.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR08.vcf > FCHR08.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR09.vcf > FCHR09.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR10.vcf > FCHR10.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR11.vcf > FCHR11.ps
awk 'BEGIN {FS="\t"; OFS="\t"} !/^#/ {print $3}' FCHR12.vcf > FCHR12.ps

awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR01.ps LRnG.STI.Chr01.ps > mod.LRnG.STI.Chr01.ps
```
Now let's add SNP.ID to the first collumn of all .ps output files and lets modify their name into mod.~

```bash
for i in LRnG.STI LRno.GRC LRno.GRS MR.GRC MR.GRS MRG.STI aLRG.STI aLRL.GRC aLRL.GRS; do 
    for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps $i.Chr$c.ps > mod.$i.Chr$c.ps 
    done 
done
```


Knowing the pseudo-heritability is also great. So let's do this. The pseudo-heritability is stored in the 6th line of the reml file that is generated by GWAS script

```bash
DIRECTORY="/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo/GWAS_RSA"

# Output file to store the 6th lines
OUTPUT_FILE="/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo/GWAS_RSA/heritability.txt"

# Check if output file already exists; if so, remove it to start fresh
if [ -f "$OUTPUT_FILE" ]; then
    rm "$OUTPUT_FILE"
fi
for file in "$DIRECTORY"/*.reml; do
    # Extract the 6th line and append it to the output file
    sed -n '6p' "$file" >> "$OUTPUT_FILE"
    # Optionally, append a note about which file the line was extracted from
    echo "Extracted from $(basename "$file")" >> "$OUTPUT_FILE"
done

echo "Extraction complete."
```

It is also good to have allele frequencies within the GWAS output file - so let's calculate them first per chromosome:

```bash
# Calculating Allele Frequency
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' CHR01.vcf > allele_counts_CHR01.txt
# Add Chr_Post into allele count file
awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2, $0}' allele_counts_CHR01.txt > AC_CHR01.txt
awk -F'\t' 'BEGIN{OFS="\t"} {print $0, $1, $2, $3, $4, ($6/$7)}' AC_CHR01.txt > AF_CHR01.txt

for c in $(seq -w 1 12); do
    awk 'BEGIN{FS=OFS="\t"} {print $1"_"$2, $0}' allele_counts_CHR$c.txt > AC_CHR$c.txt
    awk -F'\t' 'BEGIN{OFS="\t"} {print $0, $1, $2, $3, $4, ($6/$7)}' AC_CHR$c.txt > AF_CHR$c.txt
done

# Add allele frequency to all .ps files with the GWAS results
paste AF_CHR01.txt mod.LRnG.STI.Chr01.ps > AF_aLRG.STI.Chr01.txt

for j in LRnG.STI LRno.GRS MR.GRS aLRG.STI aLRL.GRS LRno.GRC MR.GRC MRG.STI aLRL.GRC;do
    for c in $(seq -w 1 12);do
        paste AF_CHR$c.txt mod.$j.Chr$c.ps > AF_$j.Chr$c.txt
    done
done
```

Now let's combine all of the .ps output from GWAS containing the p-values, effect size and so on into one file for all chromosomes:
```bash
for j in LRnG.STI LRno.GRS MR.GRS aLRG.STI aLRL.GRS LRno.GRC MR.GRC MRG.STI aLRL.GRC;do
    cat AF_$j.Chr01.txt AF_$j.Chr02.txt AF_$j.Chr03.txt AF_$j.Chr04.txt AF_$j.Chr05.txt AF_$j.Chr06.txt AF_$j.Chr07.txt AF_$j.Chr08.txt AF_$j.Chr09.txt AF_$j.Chr10.txt AF_$j.Chr11.txt AF_$j.Chr12.txt > FINAL_AF_$j.txt
done
```

Now let's re-do the remainig traits:
```bash
for e in CoG.3C LRno.4S TRS.4C CoG.3S MRL.3C TRS.4S CoG.4C MRL.3S aLRL.3C CoG.4S MRL.4C aLRL.3S LRL.3C MRL.4S trait.aLRL.4C  LRL.3S MRLpTRS.3C aLRL.4S LRL.4C MRLpTRS.3S aLRLpTRS.3C LRL.4S MRLpTRS.4C aLRLpTRS.3S LRno.3C MRLpTRS.4S aLRLpTRS.4C LRno.3S TRS.3C aLRLpTRS.4S LRno.4C TRS.3S; do
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.$e.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o $e.Chr$i &
done
done

# GWAS didnt finish yesterday, let's do remainig traits:
for i in aLRL.4C aLRLpTRS.3S LRno.3C MRLpTRS.4S aLRLpTRS.4C LRno.3S TRS.3C aLRLpTRS.4S LRno.4C TRS.3S; do
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.$e.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o $e.Chr$i &
done
done

# aLRL.4C
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.aLRL.4C.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o aLRL.4C.Chr$i &
done

 # aLRLpTRS.3S
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.aLRLpTRS.3S.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o aLRLpTRS.3S.Chr$i &
done

# LRno.3C
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.LRno.3C.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o LRno.3C.Chr$i &
done

# TRS.3C
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.TRS.3C.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o TRS.3C.Chr$i &
done

# aLRLpTRS.4S
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.aLRLpTRS.4S.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o aLRLpTRS.4S.Chr$i &
done

# LRno.4C
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.LRno.4C.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o LRno.4C.Chr$i &
done

# TRS.3S
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.TRS.3S.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o TRS.3S.Chr$i &
done

for i in CoG.3C LRno.4S TRS.4C CoG.3S MRL.3C TRS.4S CoG.4C MRL.3S aLRL.3C CoG.4S MRL.4C aLRL.3S LRL.3C MRL.4S aLRL.4C  LRL.3S MRLpTRS.3C aLRL.4S LRL.4C MRLpTRS.3S aLRLpTRS.3C LRL.4S MRLpTRS.4C aLRLpTRS.3S LRno.3C MRLpTRS.4S aLRLpTRS.4C LRno.3S TRS.3C aLRLpTRS.4S LRno.4C TRS.3S; do 
    for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps $i.Chr$c.ps > mod.$i.Chr$c.ps 
    done 
done


for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps aLRL.4C.Chr$c.ps > mod.aLRL.4C.Chr$c.ps 
done 

for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.MRLpTRS.4S.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o MRLpTRS.4S.Chr$i &
done

for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps MRLpTRS.4S.Chr$c.ps > mod.MRLpTRS.4S.Chr$c.ps 
done 
 
for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.aLRLpTRS.4C.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o aLRLpTRS.4C.Chr$i &
done

for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps aLRLpTRS.4C.Chr$c.ps > mod.aLRLpTRS.4C.Chr$c.ps 
done 

for i in $(seq -w 1 12); do
    ./emmax-intel64 -v -d 10 -t CHR$i -p trait.LRno.3S.txt -k ALLCHR.filt.prune.aBN.kinf -c ALLCHR.filt.prune.PCA.txt -o LRno.3S.Chr$i &
done

for c in $(seq -w 1 12); do 
        awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[NR]=$1; next} {$1=a[FNR]; print}' FCHR$c.ps LRno.3S.Chr$c.ps > mod.LRno.3S.Chr$c.ps 
done 

for j in CoG.3C LRno.4S TRS.4C CoG.3S MRL.3C TRS.4S CoG.4C MRL.3S aLRL.3C CoG.4S MRL.4C aLRL.3S LRL.3C MRL.4S aLRL.4C  LRL.3S MRLpTRS.3C aLRL.4S LRL.4C MRLpTRS.3S aLRLpTRS.3C LRL.4S MRLpTRS.4C aLRLpTRS.3S LRno.3C MRLpTRS.4S aLRLpTRS.4C LRno.3S TRS.3C aLRLpTRS.4S LRno.4C TRS.3S;do
    for c in $(seq -w 1 12);do
        paste AF_CHR$c.txt mod.$j.Chr$c.ps > AF_$j.Chr$c.txt
    done
done

for j in CoG.3C LRno.4S TRS.4C CoG.3S MRL.3C TRS.4S CoG.4C MRL.3S aLRL.3C CoG.4S MRL.4C aLRL.3S LRL.3C MRL.4S aLRL.4C  LRL.3S MRLpTRS.3C aLRL.4S LRL.4C MRLpTRS.3S aLRLpTRS.3C LRL.4S MRLpTRS.4C aLRLpTRS.3S LRno.3C MRLpTRS.4S aLRLpTRS.4C LRno.3S TRS.3C aLRLpTRS.4S LRno.4C TRS.3S;do
    cat AF_$j.Chr01.txt AF_$j.Chr02.txt AF_$j.Chr03.txt AF_$j.Chr04.txt AF_$j.Chr05.txt AF_$j.Chr06.txt AF_$j.Chr07.txt AF_$j.Chr08.txt AF_$j.Chr09.txt AF_$j.Chr10.txt AF_$j.Chr11.txt AF_$j.Chr12.txt > FINAL_AF_$j.txt
done

# EXTRACT HERITABILITY INTO A FILE
DIRECTORY="/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo"

# Output file to store the 6th lines
OUTPUT_FILE="/home/magdaj/Magdas_Pimp_SNPs/PimpsSNPS/per_chromo/GWAS_RSA/RSA_heritability2.txt"

for file in "$DIRECTORY"/*.reml; do
    # Extract the 6th line and append it to the output file
    sed -n '6p' "$file" >> "$OUTPUT_FILE"
    # Optionally, append a note about which file the line was extracted from
    echo "Extracted from $(basename "$file")" >> "$OUTPUT_FILE"
done
```

OK Now - let's calculate the estimated FDR per file:

```R
library(qqman)
library(data.table)
results <- read.table("FINAL_AF_CoG.3C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P.value, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))

jpeg("MAF5_FDR_CoG.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

jpeg("MAF5_ES_CoG.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()

results$trait <- "CoG.3C"
results$LOD <- -log10(results$P)
Results.sig <- subset(results, results$LOD > 4)
```

OK - let's do this for the remaining traits:

```R
results <- read.table("FINAL_AF_CoG.3S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_CoG.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_CoG.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "CoG.3S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)


# COG 4C
results <- read.table("FINAL_AF_CoG.4C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_CoG.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_CoG.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "CoG.4C"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# COG 4S
results <- read.table("FINAL_AF_CoG.4S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_CoG.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_CoG.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "CoG.4S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRL 3C
results <- read.table("FINAL_AF_LRL.3C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRL.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRL.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRL.3C"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRL 3S
results <- read.table("FINAL_AF_LRL.3S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRL.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRL.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRL.3S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRL 4C
results <- read.table("FINAL_AF_LRL.4C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRL.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRL.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRL.4C"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRL 4S
results <- read.table("FINAL_AF_LRL.4S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRL.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRL.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRL.4S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.GR STI
results <- read.table("FINAL_AF_LRnG.STI.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.GR.STI.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.GR.STI.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.GR.STI"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.3C
results <- read.table("FINAL_AF_LRno.3C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.3C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.3C"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.3S
results <- read.table("FINAL_AF_LRno.3S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.3S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.3S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.4C
results <- read.table("FINAL_AF_LRno.4C.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.4C.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.4C"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.4S
results <- read.table("FINAL_AF_LRno.4S.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.4S.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.4S"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.GR C
results <- read.table("FINAL_AF_LRno.GRC.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.GRC.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.GRC.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.GRC"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)

# LRn.GR S
results <- read.table("FINAL_AF_LRno.GRS.txt", header=F, sep="\t", fill = T)
results  <- results[,c(8:12,14:16)]
colnames(results) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
results <- na.omit(results)
results <- subset(results, results$MAF > 0.05)
results$FDR <- p.adjust(results$P, method = "BH")
results$BETA <- as.numeric(as.character(results$BETA))
jpeg("MAF5_FDR_LRno.GRS.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "FDR", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
jpeg("MAF5_ES_LRno.GRS.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(results, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "BETA", logp = TRUE,cex=0.3, ylim=c(-1,6), 
          cex.axis = 0.5,cex.main = 1)
dev.off()
results$trait <- "LRno.GRS"
results$LOD <- -log10(results$P)
Temp <- subset(results, results$LOD > 4)
Results.sig <- rbind(Results.sig, Temp)
```


## Follow up on most interesting loci

Let's explore the most interesting loci:

### LRno.GRS.CHR12

```R
library(qqman)
library(data.table)
LRno.GRS <- read.table("FINAL_AF_LRno.GRS.txt", header=F, sep="\t", fill = T)
LRno.GRS  <- LRno.GRS[,c(8:12,14:16)]
colnames(LRno.GRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRno.GRS <- na.omit(LRno.GRS)
LRno.GRS <- subset(LRno.GRS, LRno.GRS$MAF > 0.05)
head(LRno.GRS)
LRno.GRS12 <- subset(LRno.GRS, LRno.GRS$CHR == 12)
LRno.GRS12$LOD <- -log10(LRno.GRS12$P)
LRno.GRS12 <- subset(LRno.GRS12, LRno.GRS12$LOD > 5)
dim(LRno.GRS12)
min(LRno.GRS12$BP)
max(LRno.GRS12$BP) 

LRno.GRS <- read.table("FINAL_AF_LRno.GRS.txt", header=F, sep="\t", fill = T)
LRno.GRS  <- LRno.GRS[,c(8:12,14:16)]
colnames(LRno.GRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRno.GRS <- na.omit(LRno.GRS)
LRno.GRS12 <- subset(LRno.GRS, LRno.GRS$CHR == 12)
jpeg("Locus_LRnoGRS.MAF5_CHR12.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(LRno.GRS12, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
```

### MRL.GRS.CHR11

```R
MRL.GRS <- read.table("FINAL_AF_MR.GRS.txt", header=F, sep="\t", fill = T)
MRL.GRS  <- MRL.GRS[,c(8:12,14:16)]
colnames(MRL.GRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL.GRS <- na.omit(MRL.GRS)
MRL.GRS <- subset(MRL.GRS, MRL.GRS$MAF > 0.05)
head(MRL.GRS)
MRL.GRS11 <- subset(MRL.GRS, MRL.GRS$CHR == 11)
MRL.GRS11$LOD <- -log10(MRL.GRS11$P)
MRL.GRS11 <- subset(MRL.GRS11, MRL.GRS11$LOD > 5)
dim(MRL.GRS11)
MRL.GRS11
min(MRL.GRS11$BP)
max(MRL.GRS11$BP) 

MRL.GRS <- read.table("FINAL_AF_MR.GRS.txt", header=F, sep="\t", fill = T)
MRL.GRS  <- MRL.GRS[,c(8:12,14:16)]
colnames(MRL.GRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL.GRS <- na.omit(MRL.GRS)
MRL.GRS11 <- subset(MRL.GRS, MRL.GRS$CHR == 11)
jpeg("Locus_MRGRS.MAF5_CHR11.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(MRL.GRS11, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
```
### MRL.d3C CHR11

```R
MRLd3C <- read.table("FINAL_AF_MRL.3C.txt", header=F, sep="\t", fill = T)
MRLd3C  <- MRLd3C[,c(8:12,14:16)]
colnames(MRLd3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLd3C <- na.omit(MRLd3C)
MRLd3C11 <- subset(MRLd3C, MRLd3C$CHR == 11)
jpeg("Locus_MRLd3C.MAF5_CHR11.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(MRLd3C11, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()
MRLd3C11 <- subset(MRLd3C, MRLd3C$CHR == 11)
MRLd3C11$LOD <- -log10(MRLd3C11$P)
#MRLd3C11 <- subset(MRLd3C11, MRLd3C11$MAF > 0.05)
MRLd3C11 <- subset(MRLd3C11, MRLd3C11$LOD > 4)
dim(MRLd3C11)
MRLd3C11
min(MRLd3C11$BP)
max(MRLd3C11$BP) 
```

### MRLpTRS.d3S CHR9
```R
MRLpTRSd3S <- read.table("FINAL_AF_MRLpTRS.3S.txt", header=F, sep="\t", fill = T)
MRLpTRSd3S  <- MRLpTRSd3S[,c(8:12,14:16)]
colnames(MRLpTRSd3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRSd3S <- na.omit(MRLpTRSd3S)
MRLpTRSd3SC9 <- subset(MRLpTRSd3S, MRLpTRSd3S$CHR == 9)
jpeg("Locus_MRLpTRSd3S.MAF5_CHR9.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(MRLpTRSd3SC9, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#MRLpTRSd3SC9 <- subset(MRLpTRSd3SC9, MRLpTRSd3SC9$MAF > 0.05)
MRLpTRSd3SC9$LOD <- -log10(MRLpTRSd3SC9$P)
MRLpTRSd3SC9 <- subset(MRLpTRSd3SC9, MRLpTRSd3SC9$LOD > 4)
dim(MRLpTRSd3SC9)
MRLpTRSd3SC9
```

### MRLpTRS.d4S CHR9
```R
MRLpTRSd4S <- read.table("FINAL_AF_MRLpTRS.4S.txt", header=F, sep="\t", fill = T)
MRLpTRSd4S  <- MRLpTRSd4S[,c(8:12,14:16)]
colnames(MRLpTRSd4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRSd4S <- na.omit(MRLpTRSd4S)
MRLpTRSd4SC9 <- subset(MRLpTRSd4S, MRLpTRSd4S$CHR == 9)
jpeg("Locus_MRLpTRSd4S.MAF5_CHR9.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(MRLpTRSd4SC9, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#MRLpTRSd4SC9 <- subset(MRLpTRSd4SC9, MRLpTRSd4SC9$MAF > 0.05)
MRLpTRSd4SC9$LOD <- -log10(MRLpTRSd4SC9$P)
MRLpTRSd4SC9 <- subset(MRLpTRSd4SC9, MRLpTRSd4SC9$LOD > 4)
dim(MRLpTRSd4SC9)
MRLpTRSd4SC9
```

### LRno.d3S CHR2
```R
LRnod3S <- read.table("FINAL_AF_LRno.3S.txt", header=F, sep="\t", fill = T)
LRnod3S  <- LRnod3S[,c(8:12,14:16)]
colnames(LRnod3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnod3S <- na.omit(LRnod3S)
LRnod3SC2 <- subset(LRnod3S, LRnod3S$CHR == 2)
jpeg("Locus_LRnod3S.MAF5_CHR2.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(LRnod3SC2, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

LRnod3SC2 <- subset(LRnod3SC2, LRnod3SC2$MAF > 0.05)
LRnod3SC2$LOD <- -log10(LRnod3SC2$P)
LRnod3SC2 <- subset(LRnod3SC2, LRnod3SC2$LOD > 5)
dim(LRnod3SC2)
LRnod3SC2
```
### LRno.d4S CHR2
```R
LRnod4S <- read.table("FINAL_AF_LRno.4S.txt", header=F, sep="\t", fill = T)
LRnod4S  <- LRnod4S[,c(8:12,14:16)]
colnames(LRnod4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnod4S <- na.omit(LRnod4S)
LRnod4SC2 <- subset(LRnod4S, LRnod4S$CHR == 2)
jpeg("Locus_LRnod4S.MAF5_CHR2.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(LRnod4SC2, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$MAF > 0.05)
LRnod4SC2$LOD <- -log10(LRnod4SC2$P)
LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$LOD > 5)
dim(LRnod4SC2)
LRnod4SC2
```

### TRS.d3S.CHR2
```R
TRSd3S <- read.table("FINAL_AF_TRS.3S.txt", header=F, sep="\t", fill = T)
TRSd3S  <- TRSd3S[,c(8:12,14:16)]
colnames(TRSd3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRSd3S <- na.omit(TRSd3S)
TRSd3SC2 <- subset(TRSd3S, TRSd3S$CHR == 2)
jpeg("Locus_TRSd3S.MAF5_CHR2.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(TRSd3SC2, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$MAF > 0.05)
TRSd3SC2$LOD <- -log10(TRSd3SC2$P)
TRSd3SC2 <- subset(TRSd3SC2, TRSd3SC2$LOD > 4)
dim(TRSd3SC2)
TRSd3SC2
```

### TRS.d3S.CHR5
```R
TRSd3SC5 <- subset(TRSd3S, TRSd3S$CHR == 5)
jpeg("Locus_TRSd3S.MAF5_CHR5.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(TRSd3SC5, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$MAF > 0.05)
TRSd3SC5$LOD <- -log10(TRSd3SC5$P)
TRSd3SC5 <- subset(TRSd3SC5, TRSd3SC5$LOD > 4)
dim(TRSd3SC5)
TRSd3SC5
```

### TRS.d4S.CHR2
```R
TRSd4S <- read.table("FINAL_AF_TRS.4S.txt", header=F, sep="\t", fill = T)
TRSd4S  <- TRSd4S[,c(8:12,14:16)]
colnames(TRSd4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRSd4S <- na.omit(TRSd4S)
TRSd4SC2 <- subset(TRSd4S, TRSd4S$CHR == 2)
jpeg("Locus_TRSd4S.MAF5_CHR2.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(TRSd4SC2, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$MAF > 0.05)
TRSd4SC2$LOD <- -log10(TRSd4SC2$P)
TRSd4SC2 <- subset(TRSd4SC2, TRSd4SC2$LOD > 4)
dim(TRSd4SC2)
TRSd4SC2
```

### TRS.d4S.CHR5
```R
TRSd4SC5 <- subset(TRSd4S, TRSd4S$CHR == 5)
jpeg("Locus_TRSd4S.MAF5_CHR5.GWAS.jpeg", height = 8, width = 16, units="cm", res = 600)
manhattan(TRSd4SC5, chr = "CHR", bp = "BP",
          snp = "SNP",col = c("blue4", "orange3"),
          p = "P", logp = TRUE,cex=0.3, 
          cex.axis = 0.5,cex.main = 1)
dev.off()

#LRnod4SC2 <- subset(LRnod4SC2, LRnod4SC2$MAF > 0.05)
TRSd4SC5$LOD <- -log10(TRSd4SC5$P)
TRSd4SC5 <- subset(TRSd4SC5, TRSd4SC5$LOD > 4)
dim(TRSd4SC5)
TRSd4SC5
```

OK - but what would happen if I select all of the files that we have and examine overlapping associations above LOD4? 

```R
LRL3C <- read.table("FINAL_AF_LRL.3C.txt", header=F, sep="\t", fill = T)
LRL3C  <- LRL3C[,c(8:12,14:16)]
colnames(LRL3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRL3C <- na.omit(LRL3C)
LRL3C$LOD <- -log10(LRL3C$P)
LRL3C <- subset(LRL3C, LRL3C$LOD > 4)
LRL3C$trait <- "LRL3C"

LRL3S <- read.table("FINAL_AF_LRL.3S.txt", header=F, sep="\t", fill = T)
LRL3S  <- LRL3S[,c(8:12,14:16)]
colnames(LRL3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRL3S <- na.omit(LRL3S)
LRL3S$LOD <- -log10(LRL3S$P)
LRL3S <- subset(LRL3S, LRL3S$LOD > 4)
LRL3S$trait <- "LRL3S"

LRL4C <- read.table("FINAL_AF_LRL.4C.txt", header=F, sep="\t", fill = T)
LRL4C  <- LRL4C[,c(8:12,14:16)]
colnames(LRL4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRL4C <- na.omit(LRL4C)
LRL4C$LOD <- -log10(LRL4C$P)
LRL4C <- subset(LRL4C, LRL4C$LOD > 4)
LRL4C$trait <- "LRL4C"

LRL4S <- read.table("FINAL_AF_LRL.4S.txt", header=F, sep="\t", fill = T)
LRL4S  <- LRL4S[,c(8:12,14:16)]
colnames(LRL4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRL4S <- na.omit(LRL4S)
LRL4S$LOD <- -log10(LRL4S$P)
LRL4S <- subset(LRL4C, LRL4S$LOD > 4)
LRL4S$trait <- "LRL4S"

RSA <- rbind(LRL3C, LRL4C, LRL3S, LRL4S)
write.csv(RSA, "GWAS_all_LOD4.csv", row.names = F)
```
let's add more files to it:

```R
LRnoSTI <- read.table("FINAL_AF_LRnG.STI.txt" , header=F, sep="\t", fill = T)
LRnoSTI  <- LRnoSTI[,c(8:12,14:16)]
colnames(LRnoSTI) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnoSTI <- na.omit(LRnoSTI)
LRnoSTI$LOD <- -log10(LRnoSTI$P)
LRnoSTI <- subset(LRnoSTI, LRnoSTI$LOD > 4)
LRnoSTI$trait <- "LRnoSTI"
RSA <- rbind(RSA, LRnoSTI)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnon3C <- read.table("FINAL_AF_LRno.3C.txt", header=F, sep="\t", fill = T)
LRnon3C  <- LRnon3C[,c(8:12,14:16)]
colnames(LRnon3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnon3C <- na.omit(LRnon3C)
LRnon3C$LOD <- -log10(LRnon3C$P)
LRnon3C <- subset(LRnon3C, LRnon3C$LOD > 4)
LRnon3C$trait <- "LRnon3C"
RSA <- rbind(RSA, LRnon3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnon3S <- read.table("FINAL_AF_LRno.3S.txt", header=F, sep="\t", fill = T)
LRnon3S  <- LRnon3S[,c(8:12,14:16)]
colnames(LRnon3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnon3S <- na.omit(LRnon3S)
LRnon3S$LOD <- -log10(LRnon3S$P)
LRnon3S <- subset(LRnon3S, LRnon3S$LOD > 4)
LRnon3S$trait <- "LRnon3S"
RSA <- rbind(RSA, LRnon3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnon4C <- read.table("FINAL_AF_LRno.4C.txt", header=F, sep="\t", fill = T)
LRnon4C  <- LRnon4C[,c(8:12,14:16)]
colnames(LRnon4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnon4C <- na.omit(LRnon4C)
LRnon4C$LOD <- -log10(LRnon4C$P)
LRnon4C <- subset(LRnon4C, LRnon4C$LOD > 4)
LRnon4C$trait <- "LRnon4C"
RSA <- rbind(RSA, LRnon4C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnon4S <- read.table("FINAL_AF_LRno.4S.txt", header=F, sep="\t", fill = T)
LRnon4S  <- LRnon4S[,c(8:12,14:16)]
colnames(LRnon4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnon4S <- na.omit(LRnon4S)
LRnon4S$LOD <- -log10(LRnon4S$P)
LRnon4S <- subset(LRnon4S, LRnon4S$LOD > 4)
LRnon4S$trait <- "LRnon4S"
RSA <- rbind(RSA, LRnon4S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnonGRC <- read.table("FINAL_AF_LRno.GRC.txt", header=F, sep="\t", fill = T)
LRnonGRC  <- LRnonGRC[,c(8:12,14:16)]
colnames(LRnonGRC) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnonGRC <- na.omit(LRnonGRC)
LRnonGRC$LOD <- -log10(LRnonGRC$P)
LRnonGRC <- subset(LRnonGRC, LRnonGRC$LOD > 4)
LRnonGRC$trait <- "LRnonGRC"
RSA <- rbind(RSA, LRnonGRC)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

LRnonGRS <- read.table("FINAL_AF_LRno.GRS.txt", header=F, sep="\t", fill = T)
LRnonGRS  <- LRnonGRS[,c(8:12,14:16)]
colnames(LRnonGRSLRnonGRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
LRnonGRS <- na.omit(LRnonGRS)
LRnonGRS$LOD <- -log10(LRnonGRS$P)
LRnonGRS <- subset(LRnonGRS, LRnonGRS$LOD > 4)
LRnonGRS$trait <- "LRnonGRS"
RSA <- rbind(RSA, LRnonGRS)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLGRC <- read.table("FINAL_AF_MR.GRC.txt", header=F, sep="\t", fill = T)
MRLGRC  <- MRLGRC[,c(8:12,14:16)]
colnames(MRLGRC) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLGRC <- na.omit(MRLGRC)
MRLGRC$LOD <- -log10(MRLGRC$P)
MRLGRC <- subset(MRLGRC, MRLGRC$LOD > 4)
MRLGRC$trait <- "MRLGRC"
RSA <- rbind(RSA, MRLGRC)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLGRS <- read.table("FINAL_AF_MR.GRS.txt", header=F, sep="\t", fill = T)
MRLGRS  <- MRLGRS[,c(8:12,14:16)]
colnames(MRLGRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLGRS <- na.omit(MRLGRS)
MRLGRS$LOD <- -log10(MRLGRS$P)
MRLGRS <- subset(MRLGRS, MRLGRS$LOD > 4)
MRLGRS$trait <- "MRLGRS"
RSA <- rbind(RSA, MRLGRS)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLSTI <- read.table("FINAL_AF_MRG.STI.txt", header=F, sep="\t", fill = T)
MRLSTI  <- MRLSTI[,c(8:12,14:16)]
colnames(MRLSTI) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLSTI <- na.omit(MRLSTI)
MRLSTI$LOD <- -log10(MRLSTI$P)
MRLSTI <- subset(MRLSTI, MRLSTI$LOD > 4)
MRLSTI$trait <- "MRLSTI"
RSA <- rbind(RSA, MRLSTI)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)
```

And more:

```R
MRL3C <- read.table("FINAL_AF_MRL.3C.txt", header=F, sep="\t", fill = T)
MRL3C  <- MRL3C[,c(8:12,14:16)]
colnames(MRL3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL3C <- na.omit(MRL3C)
MRL3C$LOD <- -log10(MRL3C$P)
MRL3C <- subset(MRL3C, MRL3C$LOD > 4)
MRL3C$trait <- "MRL3C"
RSA <- rbind(RSA, MRL3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRL3S <- read.table("FINAL_AF_MRL.3S.txt", header=F, sep="\t", fill = T)
MRL3S  <- MRL3S[,c(8:12,14:16)]
colnames(MRL3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL3S <- na.omit(MRL3S)
MRL3S$LOD <- -log10(MRL3S$P)
MRL3S <- subset(MRL3S, MRL3S$LOD > 4)
MRL3S$trait <- "MRL3S"
RSA <- rbind(RSA, MRL3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRL4C <- read.table("FINAL_AF_MRL.4C.txt", header=F, sep="\t", fill = T)
MRL4C  <- MRL4C[,c(8:12,14:16)]
colnames(MRL4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL4C <- na.omit(MRL4C)
MRL4C$LOD <- -log10(MRL4C$P)
MRL4C <- subset(MRL4C, MRL4C$LOD > 4)
MRL4C$trait <- "MRL4C"
RSA <- rbind(RSA, MRL4C)
write.csv(RSA, "GWAS_all_LOD4.csv", row.names = F)

MRL4S <- read.table("FINAL_AF_MRL.4S.txt", header=F, sep="\t", fill = T)
MRL4S  <- MRL4S[,c(8:12,14:16)]
colnames(MRL4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRL4S <- na.omit(MRL4S)
MRL4S$LOD <- -log10(MRL4S$P)
MRL4S <- subset(MRL4S, MRL4S$LOD > 4)
MRL4S$trait <- "MRL4S"
RSA <- rbind(RSA, MRL4S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLpTRS3C <- read.table("FINAL_AF_MRLpTRS.3C.txt", header=F, sep="\t", fill = T)
MRLpTRS3C  <- MRLpTRS3C[,c(8:12,14:16)]
colnames(MRLpTRS3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRS3C <- na.omit(MRLpTRS3C)
MRLpTRS3C$LOD <- -log10(MRLpTRS3C$P)
MRLpTRS3C <- subset(MRLpTRS3C, MRLpTRS3C$LOD > 4)
MRLpTRS3C$trait <- "MRLpTRS3C"
RSA <- rbind(RSA, MRLpTRS3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLpTRS3S <- read.table("FINAL_AF_MRLpTRS.3S.txt", header=F, sep="\t", fill = T)
MRLpTRS3S  <- MRLpTRS3S[,c(8:12,14:16)]
colnames(MRLpTRS3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRS3S <- na.omit(MRLpTRS3S)
MRLpTRS3S$LOD <- -log10(MRLpTRS3S$P)
MRLpTRS3S <- subset(MRLpTRS3S, MRLpTRS3S$LOD > 4)
MRLpTRS3S$trait <- "MRLpTRS3S"
RSA <- rbind(RSA, MRLpTRS3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLpTRS4C <- read.table("FINAL_AF_MRLpTRS.4C.txt", header=F, sep="\t", fill = T)
MRLpTRS4C  <- MRLpTRS4C[,c(8:12,14:16)]
colnames(MRLpTRS4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRS4C <- na.omit(MRLpTRS4C)
MRLpTRS4C$LOD <- -log10(MRLpTRS4C$P)
MRLpTRS4C <- subset(MRLpTRS4C, MRLpTRS4C$LOD > 4)
MRLpTRS4C$trait <- "MRLpTRS4C"
RSA <- rbind(RSA, MRLpTRS4C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

MRLpTRS4S <- read.table("FINAL_AF_MRLpTRS.4S.txt", header=F, sep="\t", fill = T)
MRLpTRS4S  <- MRLpTRS4S[,c(8:12,14:16)]
colnames(MRLpTRS4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
MRLpTRS4S <- na.omit(MRLpTRS4S)
MRLpTRS4S$LOD <- -log10(MRLpTRS4S$P)
MRLpTRS4S <- subset(MRLpTRS4S, MRLpTRS4S$LOD > 4)
MRLpTRS4S$trait <- "MRLpTRS4S"
RSA <- rbind(RSA, MRLpTRS4S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

TRS3C <- read.table("FINAL_AF_TRS.3C.txt", header=F, sep="\t", fill = T)
TRS3C  <- TRS3C[,c(8:12,14:16)]
colnames(TRS3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRS3C <- na.omit(TRS3C)
TRS3C$LOD <- -log10(TRS3C$P)
TRS3C <- subset(TRS3C, TRS3C$LOD > 4)
TRS3C$trait <- "TRS3C"
RSA <- rbind(RSA, TRS3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

TRS3S <- read.table("FINAL_AF_TRS.3S.txt", header=F, sep="\t", fill = T)
TRS3S  <- TRS3S[,c(8:12,14:16)]
colnames(TRS3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRS3S <- na.omit(TRS3S)
TRS3S$LOD <- -log10(TRS3S$P)
TRS3S <- subset(TRS3S, TRS3S$LOD > 4)
TRS3S$trait <- "TRS3S"
RSA <- rbind(RSA, TRS3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

TRS4C <- read.table("FINAL_AF_TRS.4C.txt", header=F, sep="\t", fill = T)
TRS4C  <- TRS4C[,c(8:12,14:16)]
colnames(TRS4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRS4C <- na.omit(TRS4C)
TRS4C$LOD <- -log10(TRS4C$P)
TRS4C <- subset(TRS4C, TRS4C$LOD > 4)
TRS4C$trait <- "TRS4C"
RSA <- rbind(RSA, TRS4C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

TRS4S <- read.table("FINAL_AF_TRS.4S.txt", header=F, sep="\t", fill = T)
TRS4S  <- TRS4S[,c(8:12,14:16)]
colnames(TRS4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
TRS4S <- na.omit(TRS4S)
TRS4S$LOD <- -log10(TRS4S$P)
TRS4S <- subset(TRS4S, TRS4S$LOD > 4)
TRS4S$trait <- "TRS4S"
RSA <- rbind(RSA, TRS4S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLSTI <- read.table("FINAL_AF_aLRG.STI.txt", header=F, sep="\t", fill = T)
aLRLSTI  <- aLRLSTI[,c(8:12,14:16)]
colnames(aLRLSTI) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLSTI <- na.omit(aLRLSTI)
aLRLSTI$LOD <- -log10(aLRLSTI$P)
aLRLSTI <- subset(aLRLSTI, aLRLSTI$LOD > 4)
aLRLSTI$trait <- "aLRLSTI"
RSA <- rbind(RSA, aLRLSTI)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRL3C <- read.table("FINAL_AF_aLRL.3C.txt", header=F, sep="\t", fill = T)
aLRL3C  <- aLRL3C[,c(8:12,14:16)]
colnames(aLRL3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRL3C <- na.omit(aLRL3C)
aLRL3C$LOD <- -log10(aLRL3C$P)
aLRL3C <- subset(aLRL3C, aLRL3C$LOD > 4)
aLRL3C$trait <- "aLRL3C"
RSA <- rbind(RSA, aLRL3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRL3S <- read.table("FINAL_AF_aLRL.3S.txt", header=F, sep="\t", fill = T)
aLRL3S  <- aLRL3S[,c(8:12,14:16)]
colnames(aLRL3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRL3S <- na.omit(aLRL3S)
aLRL3S$LOD <- -log10(aLRL3S$P)
aLRL3S <- subset(aLRL3S, aLRL3S$LOD > 4)
aLRL3S$trait <- "aLRL3S"
RSA <- rbind(RSA, aLRL3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRL4C <- read.table("FINAL_AF_aLRL.4C.txt", header=F, sep="\t", fill = T)
aLRL4C  <- aLRL4C[,c(8:12,14:16)]
colnames(aLRL4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRL4C <- na.omit(aLRL4C)
aLRL4C$LOD <- -log10(aLRL4C$P)
aLRL4C <- subset(aLRL4C, aLRL4C$LOD > 4)
aLRL4C$trait <- "aLRL4C"
RSA <- rbind(RSA, aLRL4C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRL4S <- read.table("FINAL_AF_aLRL.4S.txt", header=F, sep="\t", fill = T)
aLRL4S  <- aLRL4S[,c(8:12,14:16)]
colnames(aLRL4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRL4S <- na.omit(aLRL4S)
aLRL4S$LOD <- -log10(aLRL4S$P)
aLRL4S <- subset(aLRL4S, aLRL4S$LOD > 4)
aLRL4S$trait <- "aLRL4S"
RSA <- rbind(RSA, aLRL4S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLGRC <- read.table("FINAL_AF_aLRL.GRC.txt", header=F, sep="\t", fill = T)
aLRLGRC  <- aLRLGRC[,c(8:12,14:16)]
colnames(aLRLGRC) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLGRC <- na.omit(aLRLGRC)
aLRLGRC$LOD <- -log10(aLRLGRC$P)
aLRLGRC <- subset(aLRLGRC, aLRLGRC$LOD > 4)
aLRLGRC$trait <- "aLRLGRC"
RSA <- rbind(RSA, aLRLGRC)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLGRS <- read.table("FINAL_AF_aLRL.GRS.txt", header=F, sep="\t", fill = T)
aLRLGRS  <- aLRLGRS[,c(8:12,14:16)]
colnames(aLRLGRS) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLGRS <- na.omit(aLRLGRS)
aLRLGRS$LOD <- -log10(aLRLGRS$P)
aLRLGRS <- subset(aLRLGRS, aLRLGRC$LOD > 4)
aLRLGRS$trait <- "aLRLGRS"
RSA <- rbind(RSA, aLRLGRS)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLpTRS3C <- read.table("FINAL_AF_aLRLpTRS.3C.txt", header=F, sep="\t", fill = T)
aLRLpTRS3C  <- aLRLpTRS3C[,c(8:12,14:16)]
colnames(aLRLpTRS3C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLpTRS3C <- na.omit(aLRLpTRS3C)
aLRLpTRS3C$LOD <- -log10(aLRLpTRS3C$P)
aLRLpTRS3C <- subset(aLRLpTRS3C, aLRLGRC$LOD > 4)
aLRLpTRS3C$trait <- "aLRLpTRS3C"
RSA <- rbind(RSA, aLRLpTRS3C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLpTRS3S <- read.table("FINAL_AF_aLRLpTRS.3S.txt", header=F, sep="\t", fill = T)
aLRLpTRS3S  <- aLRLpTRS3S[,c(8:12,14:16)]
colnames(aLRLpTRS3S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLpTRS3S <- na.omit(aLRLpTRS3S)
aLRLpTRS3S$LOD <- -log10(aLRLpTRS3S$P)
aLRLpTRS3S <- subset(aLRLpTRS3S, aLRLGRC$LOD > 4)
aLRLpTRS3S$trait <- "aLRLpTRS3S"
RSA <- rbind(RSA, aLRLpTRS3S)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLpTRS4C <- read.table("FINAL_AF_aLRLpTRS.4C.txt", header=F, sep="\t", fill = T)
aLRLpTRS4C  <- aLRLpTRS4C[,c(8:12,14:16)]
colnames(aLRLpTRS4C) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLpTRS4C <- na.omit(aLRLpTRS4C)
aLRLpTRS4C$LOD <- -log10(aLRLpTRS4C$P)
aLRLpTRS4C <- subset(aLRLpTRS4C, aLRLGRC$LOD > 4)
aLRLpTRS4C$trait <- "aLRLpTRS4C"
RSA <- rbind(RSA, aLRLpTRS4C)
write.csv(RSA, "GWAS_all_LOD4.csv", col.names = F)

aLRLpTRS4S <- read.table("FINAL_AF_aLRLpTRS.4S.txt", header=F, sep="\t", fill = T)
aLRLpTRS4S  <- aLRLpTRS4S[,c(8:12,14:16)]
colnames(aLRLpTRS4S) <- c("SNP", "CHR","BP","AL","MAF", "BETA", "SE", "P")
aLRLpTRS4S <- na.omit(aLRLpTRS4S)
aLRLpTRS4S$LOD <- -log10(aLRLpTRS4S$P)
aLRLpTRS4S <- subset(aLRLpTRS4S, aLRLGRC$LOD > 4)
aLRLpTRS4S$trait <- "aLRLpTRS4S"
RSA <- rbind(RSA, aLRLpTRS4S)
write.csv(RSA, "GWAS_all_LOD4.csv", row.names = F)

RSAMAF5 <- subsett(RSA, RSA$MAF > 0.05)
write.csv(RSAMAF5, "GWAS_all_LOD4_MAF5.csv", row.names = F)
```

```R
library(dplyr)

snp_sorted <- RSA %>% arrange(CHR, BP)

snp_sorted <- snp_sorted %>%
    group_by(CHR) %>% mutate(prev_pos = lag(BP), # Position of the previous SNP
         next_pos = lead(BP), # Position of the next SNP
         dist_to_prev = BP - prev_pos, # Distance to the previous SNP
         dist_to_next = next_pos - BP) # Distance to the next SNP

# Filtering SNPs within 10,000 bp of each other
snp_filtered <- snp_sorted %>%
  filter(dist_to_prev <= 10000 | dist_to_next <= 10000)

snp_filtered_maf <- snp_filtered %>% filter(MAF >= 0.05)

snp_unique <- unique(snp_filtered_maf$SNP, snp_filtered_maf$CHR, snp_filtered_maf$MAF)
dim(snp_unique)

snp_sorted <- snp_sorted %>%
  mutate(interval = cut(BP, breaks = seq(from = min(BP), to = max(BP) + 10000, by = 10000), labels = FALSE, include.lowest = TRUE))

# Count the number of SNPs in each interval for each chromosome
snp_counts <- snp_sorted %>%
  group_by(chr, interval) %>%
  summarise(count = n(), .groups = 'drop')
```



