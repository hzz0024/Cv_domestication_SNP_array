# Replace the original array sample id with read sample id
cd /workdir/hz269/domestication_600K/00_vcf
./1_1_snp_array_format.sh

bcftools query -l genetyped_data_all_samples.vcf > sample_original_name # 842 samples
bcftools reheader -s sample_rename genetyped_data_all_samples.vcf > genetyped_data_all_samples.rename.vcf
mv genetyped_data_all_samples.rename.vcf

 # Replace the chromosome with numbers (1-10)
./1_2_sed_chr.sh

for i in genetyped_data_all_samples.rename.vcf; do
    sed -i.bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/MT/11/g' $i
done

rm genetyped_data_all_samples.rename.vcf.bak

# edit the VCF file to include SNP ID
./1_3_change_ID_array.sh
#!/bin/sh
python3 add_array.py

cat add_array.py
fname = 'genetyped_data_all_samples.rename.vcf'
outname = fname + '.out'

idx = 0
with open(fname, 'r') as f, open(outname, 'w') as w:
    for l in f:
        if l.startswith('#'):
            pass
        else:
            idx += 1
            ss = l.split()
            chrom = ss[0]
            pos = ss[1]
            ID = chrom + '_' + pos
            if ss[2].startswith('AX'):
                ss[2] = ID
                l = '\t'.join(ss)
                l += '\n'
        w.write(l)

# exclude inversions and mtDNA
./1_4_exclude_invers.sh

vcftools --vcf genetyped_data_all_samples.rename.vcf.out --exclude-bed delly.inversions.masked.bed --recode --recode-INFO-all --out genetyped_data_all_samples_nodellyinvers

# After filtering, kept 300446 out of a possible 300446 Sites

vcftools --vcf genetyped_data_all_samples_nodellyinvers.recode.vcf --exclude-bed invers.bed --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --recode --recode-INFO-all --out genetyped_data_all_samples_noinvers

# After filtering, kept 276327 out of a possible 300446 Sites

# exclude LGF, PCs, VC familes, and CBW populations
1_5_keep_population_n_539.sh

vcftools --vcf genetyped_data_all_samples_noinvers.recode.vcf --keep sample_id_539.txt --recode --recode-INFO-all --out genetyped_data_n_509

# VCFtools - 0.1.17
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf genetyped_data_all_samples_noinvers.recode.vcf
# 	--keep sample_id_539.txt
# 	--recode-INFO-all
# 	--out genetyped_data_n_539
# 	--recode

# Keeping individuals in 'keep' list
# After filtering, kept 539 out of 842 Individuals
# Outputting VCF file...
# After filtering, kept 276327 out of a possible 276327 Sites
# Run Time = 77.00 seconds

# check the missing ind
./1_6_indmiss.sh
source /programs/miniconda3/bin/activate dDocent-2.8.13
filter_missing_ind.sh genetyped_data_n_539.recode.vcf genetyped_data_n_539_indmiss

# VCFtools - 0.1.16
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf genetyped_data_n_539.recode.vcf
# 	--missing-indv
# 	--out genetyped_data_n_539_indmiss

# After filtering, kept 539 out of 539 Individuals
# Outputting Individual Missingness
# After filtering, kept 276327 out of a possible 276327 Sites
# Run Time = 15.00 seconds



#                                           Histogram of % missing data per individual
#       400 +---------------------------------------------------------------------------------------------------------+
#           |            *             +            *            +            +            +             +            |
#           |            *                          *           'totalmissing' using (bin($1,binwidth)):(1.0) ******* |
#       350 |-+          *                          *                                                               +-|
#           |            *                          *                                                                 |
#           |            *                          *                                                                 |
#       300 |-+          *                          *                                                               +-|
#           |            *                          *                                                                 |
#           |            *                          *                                                                 |
#       250 |-+          *                          *                                                               +-|
#           |            *                          *                                                                 |
#       200 |-+          *                          *                                                               +-|
#           |            *                          *                                                                 |
#           |            *                          *                                                                 |
#       150 |-+          *                          *                                                               +-|
#           |            *                          *                                                                 |
#           |            *                          *                                                                 |
#       100 |-+          *                          *                                                               +-|
#           |            *                          ***************************                                       |
#           |            *                          *                         *                                       |
#        50 |-+          *                          *                         *                                     +-|
#           |*************                          *                         *                                       |
#           |            *             +            *            +            ****************************            |
#         0 +---------------------------------------------------------------------------------------------------------+
#         0.025         0.03         0.035         0.04        0.045         0.05        0.055          0.06        0.065
#                                                        % of missing data

# The 85% cutoff would be 0.042656
# Would you like to set a different cutoff, yes or no
# yes
# Please enter new cutoff
# 0.1
# All individuals with more than 10.0% missing data will be removed.

# VCFtools - 0.1.16
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf genetyped_data_n_539.recode.vcf
# 	--remove lowDP.indv
# 	--recode-INFO-all
# 	--out genetyped_data_n_539_indmiss
# 	--recode

# Excluding individuals in 'exclude' list
# After filtering, kept 539 out of 539 Individuals
# Outputting VCF file...
# After filtering, kept 276327 out of a possible 276327 Sites
# Run Time = 61.00 seconds

# maf and missing data filtering
./1_7_maf_rate_filter.sh
vcftools --vcf genetyped_data_n_539_indmiss.recode.vcf --maf 0.05 --max-missing 0.95 --recode --recode-INFO-all --out genetyped_data_n_539_maf05_maxmiss095

# VCFtools - 0.1.16
# (C) Adam Auton and Anthony Marcketta 2009

# Parameters as interpreted:
# 	--vcf genetyped_data_n_539_indmiss.recode.vcf
# 	--recode-INFO-all
# 	--maf 0.05
# 	--max-missing 0.95
# 	--out genetyped_data_n_539_maf05_maxmiss095
# 	--recode

# After filtering, kept 539 out of 539 Individuals
# Outputting VCF file...
# After filtering, kept 147160 out of a possible 276327 Sites
# Run Time = 50.00 seconds

./1_8_pop_missing.sh
source /programs/miniconda3/bin/activate dDocent-2.8.13
pop_missing_filter.sh genetyped_data_n_539_maf05_maxmiss095.recode.vcf popmap.txt 0.95 18 genetyped_data_n_539_maf05_maxmiss095_popmiss095

./1_9_HWE_by_pop.sh
./filter_hwe_by_pop.pl -v genetyped_data_n_539_maf05_maxmiss095_popmiss095.recode.vcf -p popmap.txt -h 0.01 -c 0.5 -o genetyped_data_n_539_maf05_maxmiss095_popmiss095_hwe

# Processing population: DBW1 (31 inds)
# Processing population: DBW2 (32 inds)
# Processing population: DBX1 (32 inds)
# Processing population: DBX2 (31 inds)
# Processing population: DBX3 (31 inds)
# Processing population: LIW1 (31 inds)
# Processing population: LIW2 (30 inds)
# Processing population: MEH2 (32 inds)
# Processing population: MEW1 (30 inds)
# Processing population: MEW2 (31 inds)
# Processing population: NCW1 (32 inds)
# Processing population: NCW2 (30 inds)
# Processing population: NEH1 (32 inds)
# Processing population: NEH2 (32 inds)
# Processing population: NYH1 (30 inds)
# Processing population: UMFS (30 inds)
# Processing population: UNC1 (20 inds)
# Processing population: UNC2 (22 inds)
# Outputting results of HWE test for filtered loci to 'filtered.hwe'
# Kept 141676 of a possible 147160 loci (filtered 5484 loci)

col_gradient <- c(  "#FA812F", "#849cc1", "#1D92BD", "#8ad5d9","#FAAF08", "#0A2C86", "#FA4032",
                    "#f9476b", "#fec155","#fddae1", "#cf7fbc",  "#e2b2d6", "#e1bb94", "#e1bb94", "#fbd0a5", "#b58383")
