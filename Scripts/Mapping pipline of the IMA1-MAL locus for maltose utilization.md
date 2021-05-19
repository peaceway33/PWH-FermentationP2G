# MALTOSE UTILIZATION

# _____________________________________________________________________________________________

## TABLE OF CONTENT

[TOC]



# _____________________________________________________________________________________________

## INFO

In collaboration with Ping-Wei Ho

We have two phenotypes:
grow VS No grow in maltose
we want to check:
- presence absence of MAL genes
- copy number of MAL genes
- eventual presence of any mutation

I have:
- two parental strains reference genomes
- raw reads of 1153 segregants
- S228C reference genome as well (e.g.: MAL regulator is mutated and not functional in S288C)
I may need to annotate the reference genomes / OR DO AN EXHAUSTIVE BLAST SEARCH



### Softwares

| Name        | v.                  | date                      | comments                     |
| ----------- | ------------------- | ------------------------- | ---------------------------- |
| bedtools    | v2.29.0-14-gea5e914 |                           | coverage depth               |
| blast       | 2.2.31+             | build Jan 7 2016 23:17:17 | sequence similarity searches |
| BWA         | 0.7.17-r1194-dirty  |                           | short reads alignments       |
| cnvnator    | v0.4.1              |                           | call CNVariants              |
| FastQC      | v0.11.8             |                           | reads QC                     |
| GATK        | 4.0.11.0            |                           | call variants                |
| samtools    | 1.9                 |                           | alignments manipulations     |
| SnpEff      | 4.3t                | build 2017-11-24 10:18    | variants annotation          |
| Trimmomatic | 0.33                |                           | short reads trim             |



# _____________________________________________________________________________________________



## IDENTIFY MAL GENE IN RM11 AND YJM975

First, we want to see if some MAL genes are missing from the RM11 YJM975 genome/proteome due to assembly/annotation quality of the publicly available data. To do so, we do a sequence similarity search (BLAST) using S288C MAL genes against the proteomes and the genome. 



### protein VS protein search

```bash
# prepare protein BLAST databases
for DB in RM11.prot.fa YJM975.prot.fa; do makeblastdb -in $DB -dbtype prot; done

# protein VS protein search
for DB in RM11.prot.fa YJM975.prot.fa;
	do blastp -db $DB \
		-query S288C.prot.MAL.fa \
		-out S288C.prot.MAL.BLAST.$(basename $DB .prot.fa) \
		-outfmt 6 -evalue 10e-45 \
		-num_threads 4;
done

# retrieve MAL proteins sequences and gff positions
for DB in RM11 YJM975;
	do cut -f 2 S288C.prot.MAL.BLAST.$DB |\
   		sort -u > $DB.prot.MAL.lst;
   	perl ~/scripts/SelectList_Fasta.pl \
   		$DB.prot.fa $DB.prot.MAL.lst > $DB.prot.MAL.fa;
done
for DB in RM11 YJM975;
	do while read line;
		do grep $line ../$DB.annot.gff |\
        	grep cds;
    done < $DB.prot.MAL.lst;
done
```





### protein VS genome search

```bash
# create genome BLAST dabatase
for DB in RM11.fa YJM975.fa; do makeblastdb -in $DB -dbtype nucl; done

# protein VS genomic sequence search
for DB in RM11.fa YJM975.fa;
	do tblastn -db $DB \
		-query S288C.prot.MAL.fa \
		-out S288C.prot.MAL.tBLASTn.$(basename $DB .fa) \
		-outfmt 6 -evalue 10e-45 \
		-num_threads 4;
done

# retrieve genomic windows with significant hits to MAL genes
for DB in RM11 YJM975 S288C;
	do cut -f 2,9,10 S288C.prot.MAL.tBLASTn.$DB |\
    	sort -k1,1 -k2,2n -k3,3n |\
        uniq > $DB.genome.pos.MAL.bed;
done

### clean overlapping positions in RM11.genome.pos.MAL.bed and YJM97.genome.pos.MAL.bed by hand
```



Assemblies of RM11 and YJM975 genomes seems not to be complete, therefore for downstream analysis it is better to refer always to S288C genome, also because more maltose related genes are identified there. **[In retrospective, this may reflect the real situation fo these two genomes.]**

Maltose utilization genes identified in are (S288C position as reference):

| Chr  | GenBank   | start   | stop    | ID      | Gene                      | RM11 | YJM975 |
| ---- | --------- | ------- | ------- | ------- | ------------------------- | ---- | ------ |
| II   | NC_001134 | 800523  | 801929  | YBR297W | MAL33                     | Yes  | Yes    |
| II   | NC_001134 | 802631  | 804475  | YBR298C | MAL31                     | Yes  | Yes    |
| II   | NC_001134 | 805351  | 807105  | YBR299W | MAL32                     | Yes  | Yes    |
| IV   | NC_001136 | 5985    | 7814    | YDL247W | MPH2                      | Yes  | Yes    |
| VII  | NC_001139 | 1067222 | 1068991 | YGR287C | IMA1                      | Yes  | Yes    |
| VII  | NC_001139 | 1070293 | 1071714 | YGR288W | MAL13                     | Yes  | No     |
| VII  | NC_001139 | 1073963 | 1075813 | YGR289C | MAL11                     | ?    | ?      |
| VII  | NC_001139 | 1076599 | 1078353 | YGR292W | MAL12                     | Yes  | No     |
| IX   | NC_001141 | 16784   | 18553   | YIL172C | IMA3                      | Yes  | Yes    |
| X    | NC_001142 | 16767   | 18536   | YJL221C | IMA4                      | Yes  | Yes    |
| X    | NC_001142 | 24341   | 26086   | YJL216C | IMA5                      | Yes  | Yes    |
| X    | NC_001142 | 738008  | 739816  | YJR160C | MPH3                      | Yes  | Yes    |
| XV   | NC_001147 | 22525   | 24294   | YOL157C | IMA2                      | Yes  | Yes    |
| XVI  | NC_001148 | 931376  | 932788  | YPR196W | Maltose_responsive_factor | Yes  | Yes    |



# _____________________________________________________________________________________________



## SEGREGANTS SAMPLES

We can download the 1,151 segregant raw reads and analyse them. We would like to look for: variants (e.g.: non-sense mutation in MAL genes), copy number variants (CNVs) interesting MAL genes that can influence the growth in maltose phenotype.

We will:

- download raw reads
- trim of low quality bases
- align to S288C, RM11, YJM975 genomes
- call and annotate variants
- identify CNVs
- correlate variants/CNVs in MAL gens with growth in maltose phentoype.



### Reads QC and trim

```bash
# trim reads
for file in *_1.fastq.gz;
	do java -jar ~/bin/trimmomatic/trimmomatic.jar PE \
		-threads 16 -phred33 \
		$file $(basename $file _1.fastq.gz)_2.fastq.gz \
		$(basename $file .fastq.gz).tr.fq.gz \
		$(basename $file .fastq.gz).un.tr.fq.gz \
		$(basename $file _1.fastq.gz)_2.tr.fq.gz \
		$(basename $file _1.fastq.gz)_2.un.tr.fq.gz \
		SLIDINGWINDOW:10:20 TRAILING:20 MINLEN:30;
done

# reads QC
for file in *.gz ; do ~/bin/FastQC/fastqc -t 16 $file; done

# orphan reads in the SRR files, I have to take care of that since Trimmomatic failed in doing so
while read line;
	do perl ../fastq-remove-orphans.pl \
		-1 "${line}"_1.tr.fq -2 "${line}"_2.tr.fq;
done < <(cut -f 2 ../MAL_all_strains.lst | sed 's/"//g')
```



Segregant "873" (reads library SRR5634436) failed.



### Variant calling

```bash
# align reads to the genomes
while read line;
	do for genome in RM11 S288C YJM975;
		do ~/bin/bwa/bwa mem -t 16 \
			-K 100000000 00_REF_genomes/"${genome}".fa \
            "${line}"_1.tr.fq.gz "${line}"_2.tr.fq.gz > \
            02_ALIGN_"${genome}"/"${line}".map."${genome}".sam;
        ~/bin/samtools-1.9/samtools view -bS -@ 16 \
        	02_ALIGN_"${genome}"/"${line}".map."${genome}".sam > \
        	02_ALIGN_"${genome}"/"${line}".map."${genome}".bam;
        ~/bin/samtools-1.9/samtools sort -@ 16 \
        	02_ALIGN_"${genome}"/"${line}".map."${genome}".bam \
        	-o 02_ALIGN_"${genome}"/"${line}".map."${genome}".sort.bam;
        rm 02_ALIGN_"${genome}"/"${line}".map."${genome}".sam \
        	02_ALIGN_"${genome}"/"${line}".map."${genome}".bam;
   done;
done < <(cut -f 2 MAL_all_strains.lst | sed 's/"//g')

# index reference genomes
for i in RM11 S288C YJM975; do ~/bin/samtools-1.9/samtools faidx $i.fa; done

# mark duplicates
for DIR in 02_*;
	do for file in $DIR/*.bam;
		do ~/bin/gatk-4.0.11.0/gatk MarkDuplicates \
			--INPUT=$DIR/$file \
			--OUTPUT=$DIR/$(basename $file .bam).md.bam \
			--METRICS_FILE=$DIR/$(basename $file .bam).md.log \
			--CREATE_INDEX=true;
	done &
done

# assign read groups
for REF in RM11 S288C YJM975;
	do for file in 02_ALIGN_"${REF}"/*.sort.md.bam;
		do ~/bin/gatk-4.0.11.0/gatk AddOrReplaceReadGroups \
			--INPUT $file \
			--OUTPUT 02_ALIGN_"${REF}"/$(basename $file .bam).r.bam \
			--RGID $(basename $file .map.$REF.sort.md.bam) \
			--RGPL $(basename $file .map.$REF.sort.md.bam) \
			--RGLB $(basename $file .map.$REF.sort.md.bam) \
			--RGPU $(basename $file .map.$REF.sort.md.bam) \
			--RGSM $(basename $file .map.$REF.sort.md.bam);
	done &
done

# run HaplotypeCaller
for REF in RM11 S288C YJM975;
	do for file in ./02_ALIGN_"${REF}"/*.bam;
		do ~/bin/gatk-4.0.11.0/gatk HaplotypeCaller \
			--input $file \
			--output 03_VARIANTS_"${REF}"/$(basename $file .bam).vcf \
			-ERC GVCF --reference 00_REF_genomes/$REF.fa \
			--base-quality-score-threshold 20 \
			--min-base-quality-score 20;
	done &
done

# combine gVCFs
for REF in RM11 S288C YJM975;
	do ~/bin/gatk-4.0.11.0/gatk CombineGVCFs \
		--reference 00_REF_genomes/$REF.fa \
		$(while read line; \
			do echo $line | \
				sed 's/^/--variant 03_VARIANTS_"\$\{REF\}"\//g' |\
                sed 's/$/.map.\$REF.sort.md.r.vcf/g' |\
                tr '\n' ' '; \
          done < <(cut -f 2 MAL_pos_strains.lst | sed 's/"//g')) \
          --output 03_VARIANTS_"${REF}"_pos/$REF.combinedgVCFs.vcf &
done

for file in *.md.r.vcf; do bgzip $file; done

## Genotype gvcf
for REF in RM11 S288C YJM975;
	do ~/bin/gatk-4.0.11.0/gatk GenotypeGVCFs \
		--reference 00_REF_genomes/$REF.fa \
		--variant 03_VARIANTS_"${REF}"/$REF.combinedgVCFs.vcf \
		--output 03_VARIANTS_"${REF}"/$REF.combinedgVCFs.gen.vcf &
done
```



### Annotate variants

```bash
# annotate variants with snpEff
for REF in RM11 S288C YJM975;
	do java -jar ~/bin/snpEff/snpEff.jar eff \
		$REF $REF.combinedgVCFs.gen.MALpos.vcf >\
        $REF.combinedgVCFs.gen.MALpos.snpEff.vcf;
    mv snpEff_genes.txt $REF.combinedgVCFs.gen.MALpos.snpEff_genes.txt;
    mv snpEff_summary.html $REF.combinedgVCFs.gen.MALpos.snpEff_summary.html;
done
```



### CNVs annotation

```bash
# run CNVnator
for REF in RM11 S288C YJM975;
	do while read line;
		do docker run -v /media/DISK2-3TB/Ping/:/data wwliao/cnvnator \
			cnvnator -root ./10_CNVs_"${REF}"_pos/out.$line.root \
			-genome ./10_CNVs_"${REF}"/00_ref_genome/"${REF}".fa \
			-tree ./02_ALIGN_"${REF}"/$line.map."${REF}".sort.md.r.bam;
			for BIN in 500 1000;
				do docker run -v /media/DISK2-3TB/Ping/:/data wwliao/cnvnator \
					cnvnator -root ./10_CNVs_"${REF}"/out.$line.root \
					-genome ./10_CNVs_"${REF}"/00_ref_genome/"${REF}".fa \
					-his $BIN \
					-d ./10_CNVs_"${REF}"_pos/;
				docker run -v /media/DISK2-3TB/Ping/:/data wwliao/cnvnator \
					cnvnator -root ./10_CNVs_"${REF}"/out.$line.root \
					-stat $BIN;
				docker run -v /media/DISK2-3TB/Ping/:/data wwliao/cnvnator \
					cnvnator -ngc \
					-root ./10_CNVs_"${REF}"/out.$line.root \
					-partition $BIN;
				docker run -v /media/DISK2-3TB/Ping/:/data wwliao/cnvnator \
					cnvnator -ngc \
					-root ./10_CNVs_"${REF}"_/out.$line.root \
					-call $BIN > \
					./10_CNVs_"${REF}"/"${REF}"."${line}".CNV_"${BIN}"bin.tab;
			done;
	done < <(cut -f 2 MAL_all_strains.lst | sed 's/"//g') &
done

# merge 500bp and 1000bp windows
for REF in RM11 S288C YJM975;
	do while read line;
		do python3.5 /media/HD_DATA/CNVnator_merger.py \
			--input_1 10_CNVs_"${REF}"/$REF.$line.CNV_500bin.tab \
			--input_2 10_CNVs_"${REF}"/$REF.$line.CNV_1000bin.tab \
			--sample $REF.$line > 10_CNVs_"${REF}"/$REF.$line.CNVmerged.tab;
	done < <(cut -f 2 MAL_all_strains.lst | sed 's/"//g') &
done

# merge output for all samples
for REF in RM11 S288C YJM975;
	do mkdir 10_CNVs_"${REF}"/02_CNVoutput_merged;
	cat 10_CNVs_"${REF}"_pos/$REF.SRR*.CNVmerged.tab |\
    	sort -k1,1 -k2,2 -k3,3n -k4,4n > 10_CNVs_"${REF}"/$REF.allMALpos_samples.tab;
    mv 10_CNVs_"${REF}"/$REF.SRR* 10_CNVs_"${REF}"/02_CNVoutput_merged;
done

# merge all results together 
python3.5 /media/HD_DATA/Ping_MaltoseGenesParser.py \
	--strainlist MAL_pos_strains.lst \
	--Rlist RM11.Maltose_genes.lst \
	--Slist S288C.Maltose_genes.lst \
	--Ylist YJM975.Maltose_genes.lst \
	--Rvar 04_snpEff_RM11_pos/RM11.combinedgVCFs.gen.MALpos.snpEff.vcf \
	--Svar 04_snpEff_S288C_pos/S288C.combinedgVCFs.gen.MALpos.snpEff.vcf \
	--Yvar 04_snpEff_YJM975_pos/YJM975.combinedgVCFs.gen.MALpos.snpEff.vcf \
	--Rcnv 10_CNVs_RM11_pos/RM11.allMALpos_samples.tab \
	--Scnv 10_CNVs_S288C_pos/S288C.allMALpos_samples.tab \
	--Ycnv 10_CNVs_YJM975_pos/YJM975.allMALpos_samples.tab > \
	Ping_MaltoseGenes_pos.results.txt
```





### Maltose genes analysis



```bash
# isolate variants for all the segregants called in each maltose gene

# For some limited number of position with no reads coverage, instead of printing out them as "./.:0,0:0:.:.:.:0,0,0", GATK prints them as "./.:0.0" or "./.:0.0.0". In order to avoid downstram errors with my scripts, I need to correct for the missing qualifiers to "./.:0,0:0:.:.:.:0,0,0"

while read line;
	do GENE=$(echo $line | cut -f 5 -d " ");
	CHR=$(echo $line | cut -f 1 -d " ");
	START=$(echo $line | cut -f 2 -d " ");
	STOP=$(echo $line | cut -f 3 -d " ");
	~/bin/gatk-4.0.11.0/gatk SelectVariants \
		--variant S288C.combinedgVCFs.gen.vcf.gz \
		--output S288C.combinedgVCFs.$GENE.gen.vcf \
		--intervals "${CHR}":"${START}"-"${STOP}"; 
	sed -i "s#./.:0,0\t#./.:0,0:0:.:.:.:0,0,0\t#g" S288C.combinedgVCFs.$GENE.gen.vcf;
	sed -i "s#./.:0,0,0\t#./.:0,0:0:.:.:.:0,0,0\t#g" S288C.combinedgVCFs.$GENE.gen.vcf;
	python3.5 /media/HD_DATA/Ping_IMA1_SNPs_simplified.py \
		--input S288C.combinedgVCFs.$GENE.gen.vcf \
		> S288C.combinedgVCFs.$GENE.gen.simpl.tab;
	while read line;
		do OLD=$(echo $line | cut -d " " -f 2 | sed 's/"//g');
		NEW=$(echo $line | cut -d " " -f 1 | sed 's/"//g');
		sed -i "s/$OLD/$NEW/g" S288C.combinedgVCFs.$GENE.gen.simpl.tab;
	done < /media/DISK2-3TB/Ping/strain_file_CORRECT20200424.txt;
done < ../S288C.Maltose_genes.lst
```



#### MAL31, MAL32, MAL33

MAL31-33 locus is situated in chromosome 2 in S288C. Among the segregants. there are two MAL31 allelic variants in, but only SNV886 (associated with RM11) contributes to the phenotype. SNV886 It is a missense mutation (G -> C; Ser328Thr) and seems not to be in the putative substrate translocation pore. There are two distinct alleles for MAL32 and MAL33 at the respective loci, but they do not correlate with the phenotype.



**S288C_MAL31-33.locus.png**

![S288C_MAL31-33.locus](./Pics_Maltose/S288C_MAL31-33.locus.png)



**MAL31_maltosegrowth.png**

![](./Pics_Maltose/MAL31_maltosegrowth.png)



**MAL32_maltosegrowth.png**

![](./Pics_Maltose/MAL32_maltosegrowth.png)



**MAL33_maltosegrowth.png**

![](./Pics_Maltose/MAL33_maltosegrowth.png)





#### MPH2

Maltose Permease Homolog (Alpha-glucoside permease): transports maltose, maltotriose, alpha-methylglucoside, and turanose. Almost identical to MPH3. The locus is situated in the subtelomeric region of Chromosome IV. The gene seems not to contribute to the phenotype. Given the weird distribution of the reads between MPH2 and MPH3 I would  say that there is only one of the two genes (combining 25% + 75% we get  100% coverage for one of the two genes). The genes have 99% sequence  identity, here explained the weird mapping and the low frequency mismatches called as variants, A in MPH2  and G in MPH3).



**S288C_MPH2.locus.png**

![](./Pics_Maltose/S288C_MPH2.locus.png)



**MPH2_maltosegrowth.png**

![](./Pics_Maltose/MPH2_maltosegrowth.png)





#### IMA1, MAL11, MAL12, MAL13

This locus is located at the subtelomeric region of chromosome VII. YJM975 segregants have a shorter allele of IMA1 and misses ~ 8kb, corresponding as well to MAL13 and MAL11 positions. IMA1, MAL13 and MAL11 contribute to the phenotype. Missense mutations in IMA1 RM11 allele (SNV4838 and SNV4829) seems not to be directly involved in the catalysis. The long complete IMA1 allele present in RM11 segregants correlates with the phenotype. RM11 and YJM975 share the same MAL12 allele. However, sequence coverage drops to 0 in multiple positions in these gene, which may be completely absent or not functional.



**S288C_IMA1_MAL11-13.locus.png**

![](./Pics_Maltose/S288C_IMA1_MAL11-13.locus.png)



**IMA1_maltosegrowth.png**

![](./Pics_Maltose/IMA1_maltosegrowth.png)



**MAL11_maltosegrowth.png**

![](./Pics_Maltose/MAL11_maltosegrowth.png)



**MAL13_maltosegrowth.png**

![](./Pics_Maltose/MAL13_maltosegrowth.png)





```bash
# Calculate IMA1 locus coverage among segregants

# create appropriate 1 kb window file to calculate coverage
cut -f 2,3 RM11.dict | tail -n +2 | sed 's/LN://g' | sed 's/SN://g' > RM11.bedchor
grep CH408044 RM11.bedchor > RM11.CH408044.bedchor
~/bin/bedtools2/bin/bedtools makewindows \
	-g RM11.CH408044.bedchor -w 1000 > RM11.CH408044.1kb_win.tab
head -n 30 RM11.CH408044.1kb_win.tab > RM11.CH408044.30kb.1kb_win.tab

# calculate the coverage on RM11 IMA1 locus 
for file in /media/DISK2-3TB/Ping/02_ALIGN_RM11_pos/*.r.bam \
	/media/DISK2-3TB/Ping/02_ALIGN_RM11/*.r.bam ChimayPop_Y*.merged.r.bam;
	do ~/bin/bedtools2/bin/bedtools coverage \
		-a /media/DISK2-3TB/Ping/00_REF_genomes/RM11.CH408044.30kb.1kb_win.tab \
		-b $file -mean > $(basename $file .map.RM11.sort.md.r.bam).CH408044.30kb.1kb_cov.bed; done

# reformat output into a data matrix
python3.5 /media/HD_DATA/Ping_IMA1_CH408044.30kb_cov.py \
	--input /media/DISK2-3TB/Ping/strain_file_CORRECT20200424.txt \
	> RM11.IMA1_10-20kb_cov.lst
```



**IMA1_allele_maltosegrowth.png**

![](./Pics_Maltose/IMA1_allele_maltosegrowth.png)





#### IMA3

IMA3 locus is situated in the subtelomeric region of Chromosome IX. RM11 and YJM975 have the same allele in these locus. The gene seems not to contribute to the phenotype.



**S288C_IMA3.locus.png**

![](./Pics_Maltose/S288C_IMA3.locus.png)





#### IMA4, IMA5

IMA4 and IMA5 locus is situated in the subtelomeric region of Chromosome X. While RM11 and YJM975 have identical IMA4 alleles and no differential variants can be called between the segregants, the different alleles of IMA5 seems not to contribute to the phenotype. SNV5983 is a missense mutation (A -> G, Ile61Val).



**S288C_IMA4-5.locus.png**

![](./Pics_Maltose/S288C_IMA4-5.locus.png)



**IMA5_maltosegrowth.png**

![](./Pics_Maltose/IMA5_maltosegrowth.png)





#### MPH3

Almost identical to MPH2. The locus is situated in the subtelomeric region of Chromosome X. The gene seems not to contribute to the phenotype. 



**S288C_MPH3.locus.png**

![](./Pics_Maltose/S288C_MPH3.locus.png)



**MPH3_maltosegrowth.png**

![](./Pics_Maltose/MPH3_maltosegrowth.png)





#### IMA2

RM11 and YJM975 share the same allele for IMA2, so this gene is not responsible for the different phenotype.



**S288C_IMA2.locus.png**

![](./Pics_Maltose/S288C_IMA2.locus.png)





#### YPR196W

Maltose responsive factor located in the subtelomeric region of chromosome XVI. RM11 and YJM975 segregants carry two distinct alleles at this locus, but the gene seems not to contribute to the phenotype.



**S288C_YPR19W.locus.png**

![](./Pics_Maltose/S288C_YPR196W.locus.png)



**YPR196W_maltosegrowth.png**

![](./Pics_Maltose/YPR196W_maltosegrowth.png)



# _____________________________________________________________________________________________



## IMA1-MAL LOCUS COVERAGE



### S288C Reference

Most of the phenotype is explained by the IMA1-MAL locus at the sub-telomeric region of chromosome 7. Let's see what happens int he segregants at this locus. First we will map the reads to S288C genome (again), then calculate the coverage over small windows (100 bp). For a fair comparison between segregants, we will calculate the relative coverage: therefore, for each position, we will divide the absolute read coverage by the average coverage of chromosome 7. This is necessary because we have high variance on the read depth sequenced. 

We will then plot the coverage information of IMA1-MAL locus for each of the 1,151 segregants as heatmap. To plot two distinct gradient colors on a heatmap (one for each parental genotype, RM11 and YJM975), we have to cheat on the ggplot system, since it can handle only one scale color gradient at the time. what we do, we create two plots, one with the Y genotype (and the R genotype transparent) and overlay it with a plot with the R genotype having  Y genotype transparent. To achieve this, we have to set an additional column reporting the alpha values (alpha = 1: full color; alpha = 0: transparent).

We will add as well a scatterplot illustrating the relationship between locus coverage and growth in maltose phenotype.



```bash
## download, QC, trim, align to S288C
while read line; do
	fasterq-dump --threads 72 --split-files --progress $line;
	java -jar ~/bin/trimmomatic/trimmomatic.jar PE \
		-threads 72 -phred33 \
		"${line}"_1.fastq "${line}"_2.fastq \
		$line.R1.tr.fq.gz $line.R1.tr.un.fq.gz \
		$line.R2.tr.fq.gz $line.R2.tr.un.fq.gz \
		SLIDINGWINDOW:10:30 TRAILING:30 MINLEN:50;
	rm *fastq *.un.*;
	bwa mem -t 72 -K 100000000 \
    	../00_REF_genomes/S288C.fa \
    	"${line}".R1.tr.fq.gz "${line}".R2.tr.fq.gz > "${line}".align.sam;
    samtools view -@ 72 -Sb "${line}".align.sam > "${line}".align.bam;
    samtools sort -@ 72 "${line}".align.bam "${line}".align.sort;
    rm *align.sam *align.bam;
done < SRR_Acc_List.txt

## coverage
bedtools makewindows -g S288C.bedchr -w 100 > S288C.100bp_win.bed
grep NC_001139 S288C.100bp_win.bed > S288C.NC_001139_.100bp_win.bed
for file in ../*.align.sort.bam; do
	NAME=$(basename $file .align.sort.bam); 
	bedtools coverage -a ~/Ping/00_REF_genomes/S288C.NC_001139_.100bp_win.bed \
		-b $file -mean > $NAME.S288C."${i}"bp_win.cov.bed;
done

# calculate chr7 coverage and normalize windows coverages
for file in *.100bp_win.cov.bed; do
	NAME=$(basename $file .S288C.100bp_win.cov.bed);
	AVG=$(awk '{ total += $4; count++ } END { print total/count }' $file);
	while read line ; do
		VALUE=$(echo $line | cut -f 4 -d ' ');
		NORM=$( echo "scale=8; $VALUE/$AVG" | bc);
		echo $line | tr ' ' '\t' | sed "s/^/$NAME\t/g" | sed "s/$/\t$NORM/g";
	done < <(bedtools intersect -wb -a S288C.IMA1.win.bed -b $file | cut -f 4-8);
done > temp1

while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d ' ');
	GEN=$(echo $line | cut -f 2 -d ' ');
	grep "${SAMPLE}" temp1 | sed "s/$/\t$GEN/g";
done < sample2genotype.tab > temp2

# split dataset and add alpha values
while read line; do
	FIRST=$(echo $line | cut -f 1-4 -d ' ' | tr ' ' '\t');
	COV=$(echo $line | cut -f 5 -d ' ');
	NORM=$(echo $line | cut -f 6 -d ' ');
	GEN=$(echo $line | cut -f 8 -d ' ');
	if [ $GEN == "Y" ]; then
		echo $line$'\t'1 | tr ' ' '\t';
		echo $FIRST$'\t'0.0000000$'\t'0$'\t'B$'\t'R$'\t'0 | tr ' ' '\t';
	elif [ $GEN == "R" ]; then
		echo $line$'\t'1 | tr ' ' '\t';
		echo $FIRST$'\t'0.0000000$'\t'0$'\t'B$'\t'Y$'\t'0 | tr ' ' '\t';
	fi;
done < temp2 | sort -k 8 > temp3

cat temp3 allsamples.100bp.cov.malgrowth.YR.tab > allsamples.100bp.cov.heatmap.YR.tab
rm temp1 temp2 temp3

# generate data for dotplot
for file in *.100bp_win.cov.bed; do
	NAME=$(basename $file .S288C.100bp_win.cov.bed);
	AVG=$(awk '{ total += $4; count++ } END { print total/count }' $file);
	while read line ; do
		VALUE=$(echo $line | cut -f 4 -d ' ');
		NORM=$( echo "scale=8; $VALUE/$AVG" | bc);
		echo $NORM;
	done < <(bedtools intersect -wb -a ../S288C.IMA1.shortallele.bed -b $file | cut -f 4-8) |\
    	awk '{ total += $1; count++ } END { print total/count }' - | sed "s/^/$NAME\t/g";
done > ../allsamples.S288C.IMA1.shortallele.YR.tab

# add genotype info as third column
while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d " ");
	GEN=$(echo $line | cut -f 2 -d " ");
	grep "${SAMPLE}" allsamples.S288C.IMA1.shortallele.YR.tab.1 | sed "s/$/\t$GEN/g";
done < sample2genotype.tab > allsamples.S288C.IMA1.shortallele.YR.tab

# swap sample name with maltose growth
echo Maltose_growth$'\t'IMA1_cov$'\t'gen > allsamples.100bp.cov.correlation.YR.tab;
while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d " ");
	PHEN=$(echo $line | cut -f 2 -d " ");
	grep "${SAMPLE}" allsamples.S288C.IMA1.shortallele.YR.tab | sed "s/$SAMPLE/$PHEN/g";
done < sample2malgrowth.tab >> allsamples.100bp.cov.correlation.YR.tab
```





**S288C.IMA1cov.heatmap.R**

```R
library("dplyr")
library("ggplot2")
library("ggpubr")
library("grid")
library("gridExtra")
library("plyr")
library("scales")
library("RColorBrewer")
library("ggnewscale")


# upload and organize
setwd("/media/DISK2-3TB//Ping/31_IMA1_cov_heatmap/02_coverage/")
cov_file = read.delim("allsamples.100bp.cov.heatmap.YR.tab", header = FALSE)
correlation_file = read.delim("allsamples.100bp.cov.correlation.YR.tab", header = TRUE)
MALgenes = read.delim("strains_MAL_gen2phe.txt", header = TRUE)

cov_file$V1 = factor(cov_file$V1,
                     levels = c("SRR5634755", "SRR5634774", "SRR5630056", "SRR5634576", "SRR5630345", "SRR5634504", "SRR5630053", "SRR5630190", "SRR5630087", "SRR5634745", "SRR5629932", "SRR5634542", "SRR5630013", "SRR5634491", "SRR5634513", "SRR5634546", "SRR5634530", "SRR5630171", "SRR5634612", "SRR5630175", "SRR5634823", "SRR5630349", "SRR5634650", "SRR5634752", "SRR5629824", "SRR5634416", "SRR5634482", "SRR5629863", "SRR5629972", "SRR5634730", "SRR5630447", "SRR5634498", "SRR5630208", "SRR5634441", "SRR5630164", "SRR5634349", "SRR5634368", "SRR5634725", "SRR5634767", "SRR5629868", "SRR5634391", "SRR5634638",
                                "SRR5629859", "SRR5630263", "SRR5630419", "SRR5634667", "SRR5634507", "SRR5634666", "SRR5630353", "SRR5634460", "SRR5629879", "SRR5629872", "SRR5630424", "SRR5634547", "SRR5634657", "SRR5630200", "SRR5634759", "SRR5634661", "SRR5630362", "SRR5634676", "SRR5629810", "SRR5630147", "SRR5630299", "SRR5629845", "SRR5630331", "SRR5630356", "SRR5630399", "SRR5634466", "SRR5630052", "SRR5630103", "SRR5630291", "SRR5634408", "SRR5630093", "SRR5629950", "SRR5634532", "SRR5629838", "SRR5630368", "SRR5630076", "SRR5634463", "SRR5630211", "SRR5630040", "SRR5629851", "SRR5634421", "SRR5634645",
                                "SRR5634686", "SRR5629801", "SRR5634586", "SRR5630121", "SRR5630034", "SRR5629835", "SRR5630351", "SRR5630023", "SRR5629830", "SRR5630032", "SRR5630444", "SRR5634431", "SRR5630438", "SRR5634596", "SRR5630320", "SRR5634364", "SRR5634459", "SRR5634540", "SRR5630048", "SRR5634786", "SRR5634558", "SRR5629805", "SRR5630192", "SRR5630074", "SRR5630148", "SRR5630105", "SRR5629891", "SRR5634362", "SRR5629920", "SRR5634438", "SRR5634406", "SRR5629864", "SRR5634422", "SRR5630425", "SRR5630028", "SRR5630326", "SRR5634499", "SRR5634706", "SRR5630170", "SRR5630179", "SRR5630233", "SRR5634810",
                                "SRR5630128", "SRR5630225", "SRR5629794", "SRR5630243", "SRR5630155", "SRR5630384", "SRR5630264", "SRR5634570", "SRR5630045", "SRR5629959", "SRR5630374", "SRR5630395", "SRR5634743", "SRR5630246", "SRR5629934", "SRR5634601", "SRR5634402", "SRR5629973", "SRR5634414", "SRR5630129", "SRR5634673", "SRR5634519", "SRR5630276", "SRR5629981", "SRR5634696", "SRR5634740", "SRR5630254", "SRR5630337", "SRR5634369", "SRR5634618", "SRR5629829", "SRR5630304", "SRR5634575", "SRR5634735", "SRR5634640", "SRR5629987", "SRR5634516", "SRR5634736", "SRR5629825", "SRR5629978", "SRR5630184", "SRR5634409",
                                "SRR5634792", "SRR5629791", "SRR5629871", "SRR5634581", "SRR5634616", "SRR5634452", "SRR5634723", "SRR5634512", "SRR5634642", "SRR5630266", "SRR5630448", "SRR5634588", "SRR5634430", "SRR5629874", "SRR5634813", "SRR5634500", "SRR5630280", "SRR5630379", "SRR5634753", "SRR5630270", "SRR5629849", "SRR5630008", "SRR5630375", "SRR5629913", "SRR5634824", "SRR5634623", "SRR5630086", "SRR5630396", "SRR5630365", "SRR5634388", "SRR5630075", "SRR5634567", "SRR5634802", "SRR5634423", "SRR5629941", "SRR5634660", "SRR5634350", "SRR5634474", "SRR5630146", "SRR5630389", "SRR5630415", "SRR5634481",
                                "SRR5634634", "SRR5630242", "SRR5630406", "SRR5629837", "SRR5630257", "SRR5630031", "SRR5630065", "SRR5630322", "SRR5634505", "SRR5630080", "SRR5634539", "SRR5634456", "SRR5630376", "SRR5634415", "SRR5634665", "SRR5630113", "SRR5634476", "SRR5630067", "SRR5629806", "SRR5629886", "SRR5630361", "SRR5630039", "SRR5634494", "SRR5629847", "SRR5629788", "SRR5634424", "SRR5630324", "SRR5630369", "SRR5630259", "SRR5630329", "SRR5634502", "SRR5630012", "SRR5634508", "SRR5634681", "SRR5634777", "SRR5634744", "SRR5630152", "SRR5634687", "SRR5630081", "SRR5634708", "SRR5634535", "SRR5634428",
                                "SRR5634348", "SRR5634716", "SRR5630255", "SRR5630385", "SRR5634822", "SRR5629903", "SRR5634464", "SRR5634403", "SRR5630068", "SRR5634399", "SRR5630109", "SRR5630408", "SRR5634677", "SRR5629968", "SRR5634798", "SRR5630106", "SRR5630073", "SRR5634515", "SRR5630290", "SRR5629971", "SRR5629974", "SRR5629858", "SRR5629979", "SRR5630104", "SRR5634603", "SRR5634492", "SRR5634707", "SRR5630413", "SRR5634384", "SRR5634746", "SRR5629841", "SRR5634779", "SRR5629797", "SRR5629834", "SRR5634347", "SRR5629820", "SRR5630046", "SRR5630295", "SRR5629840", "SRR5630241", "SRR5630205", "SRR5630265",
                                "SRR5634655", "SRR5630357", "SRR5630077", "SRR5634375", "SRR5629815", "SRR5630178", "SRR5629792", "SRR5630101", "SRR5630433", "SRR5630401", "SRR5630403", "SRR5630286", "SRR5630303", "SRR5634719", "SRR5634549", "SRR5634790", "SRR5634554", "SRR5629969", "SRR5634527", "SRR5630102", "SRR5630364", "SRR5634437", "SRR5634580", "SRR5630199", "SRR5634478", "SRR5634643", "SRR5634592", "SRR5634607", "SRR5630273", "SRR5629857", "SRR5630227", "SRR5629782", "SRR5630333", "SRR5629921", "SRR5630334", "SRR5629867", "SRR5630082", "SRR5630232", "SRR5630397", "SRR5629962", "SRR5634713", "SRR5634383",
                                "SRR5630250", "SRR5629882", "SRR5630318", "SRR5634602", "SRR5634633", "SRR5630026", "SRR5630226", "SRR5630132", "SRR5630383", "SRR5634649", "SRR5629823", "SRR5634521", "SRR5634386", "SRR5630177", "SRR5629802", "SRR5634714", "SRR5629951", "SRR5634405", "SRR5629915", "SRR5634511", "SRR5634652", "SRR5634754", "SRR5634747", "SRR5630301", "SRR5630244", "SRR5629964", "SRR5634738", "SRR5630358", "SRR5634390", "SRR5630180", "SRR5629832", "SRR5634750", "SRR5634709", "SRR5630214", "SRR5630278", "SRR5634670", "SRR5629929", "SRR5630268", "SRR5630439", "SRR5634731", "SRR5629911", "SRR5634662",
                                "SRR5634401", "SRR5634587", "SRR5630251", "SRR5630445", "SRR5629984", "SRR5634366", "SRR5630112", "SRR5630138", "SRR5634598", "SRR5630411", "SRR5634632", "SRR5634473", "SRR5634704", "SRR5630123", "SRR5634694", "SRR5630169", "SRR5630154", "SRR5634804", "SRR5630159", "SRR5634659", "SRR5629908", "SRR5630118", "SRR5630017", "SRR5634577", "SRR5629963", "SRR5629822", "SRR5629997", "SRR5634639", "SRR5630441", "SRR5630094", "SRR5630410", "SRR5634664", "SRR5630248", "SRR5630370", "SRR5634674", "SRR5634668", "SRR5630328", "SRR5634720", "SRR5634726", "SRR5630224", "SRR5629831", "SRR5634520",
                                "SRR5629783", "SRR5634626", "SRR5630186", "SRR5630089", "SRR5634454", "SRR5634648", "SRR5634465", "SRR5634372", "SRR5630409", "SRR5630137", "SRR5634757", "SRR5634772", "SRR5630041", "SRR5630124", "SRR5634352", "SRR5629875", "SRR5634742", "SRR5634784", "SRR5630059", "SRR5630209", "SRR5634449", "SRR5634410", "SRR5630422", "SRR5630258", "SRR5634543", "SRR5634773", "SRR5629839", "SRR5634453", "SRR5630236", "SRR5634443", "SRR5629796", "SRR5630191", "SRR5634710", "SRR5629976", "SRR5634551", "SRR5630204", "SRR5630058", "SRR5630267", "SRR5629860", "SRR5634712", "SRR5630157", "SRR5634404",
                                "SRR5630231", "SRR5634815", "SRR5630294", "SRR5634379", "SRR5634518", "SRR5634669", "SRR5634761", "SRR5630055", "SRR5630256", "SRR5629938", "SRR5634729", "SRR5634688", "SRR5630011", "SRR5629826", "SRR5634495", "SRR5630437", "SRR5630153", "SRR5630033", "SRR5634658", "SRR5634411", "SRR5634387", "SRR5634426", "SRR5630070", "SRR5630293", "SRR5634796", "SRR5630330", "SRR5630079", "SRR5629989", "SRR5630402", "SRR5629878", "SRR5629936", "SRR5630090", "SRR5634764", "SRR5630078", "SRR5634357", "SRR5630125", "SRR5630335", "SRR5634485", "SRR5630063", "SRR5634628", "SRR5630271", "SRR5634826",
                                "SRR5629967", "SRR5634819", "SRR5629862", "SRR5630141", "SRR5634733", "SRR5629970", "SRR5629827", "SRR5630061", "SRR5629998", "SRR5630366", "SRR5630016", "SRR5630161", "SRR5630150", "SRR5634442", "SRR5630085", "SRR5634395", "SRR5634617", "SRR5630194", "SRR5630450", "SRR5634544", "SRR5634563", "SRR5634351", "SRR5629896", "SRR5629861", "SRR5629999", "SRR5629870", "SRR5634425", "SRR5634717", "SRR5634522", "SRR5630115", "SRR5630359", "SRR5634553", "SRR5634470", "SRR5630163", "SRR5634806", "SRR5630219", "SRR5630066", "SRR5630037", "SRR5630202", "SRR5634739", "SRR5634446", "SRR5634613",
                                "SRR5630285", "SRR5630350", "SRR5630183", "SRR5634365", "SRR5630160", "SRR5634705", "SRR5630394", "SRR5630310", "SRR5630088", "SRR5629866", "SRR5634417", "SRR5634809", "SRR5634797", "SRR5634615", "SRR5634698", "SRR5629994", "SRR5630313", "SRR5630292", "SRR5634487", "SRR5634367", "SRR5629892", "SRR5634557", "SRR5630393", "SRR5629808", "SRR5634461", "SRR5630340", "SRR5629865", "SRR5634600", "SRR5630346", "SRR5630317", "SRR5634629", "SRR5629807", "SRR5629965", "SRR5629986", "SRR5630275", "SRR5634778", "SRR5629947", "SRR5629996", "SRR5630314", "SRR5629977", "SRR5629955", "SRR5634702",
                                "SRR5630309", "SRR5630440", "SRR5629894", "SRR5629884", "SRR5634548", "SRR5630381", "SRR5634545", "SRR5634419", "SRR5634793", "SRR5634541", "SRR5630151", "SRR5630064", "SRR5630230", "SRR5634765", "SRR5634671", "SRR5634489", "SRR5634560", "SRR5634477", "SRR5630047", "SRR5634497", "SRR5634534", "SRR5630083", "SRR5630127", "SRR5634801", "SRR5634398", "SRR5634825", "SRR5630174", "SRR5634469", "SRR5634816", "SRR5634644", "SRR5634611", "SRR5630238", "SRR5634814", "SRR5630261", "SRR5634631", "SRR5629990", "SRR5629924", "SRR5630418", "SRR5629786", "SRR5630050", "SRR5629933", "SRR5629926",
                                "SRR5630036", "SRR5629833", "SRR5634355", "SRR5634536", "SRR5634400", "SRR5634766", "SRR5629985", "SRR5630018", "SRR5630100", "SRR5634396", "SRR5630193", "SRR5634363", "SRR5629873", "SRR5634599", "SRR5634562", "SRR5634808", "SRR5634697", "SRR5634732", "SRR5634389", "SRR5634734", "SRR5630071", "SRR5630404", "SRR5630042", "SRR5630414", "SRR5630057", "SRR5634353", "SRR5630020", "SRR5634589", "SRR5630315", "SRR5629905", "SRR5634486", "SRR5634371", "SRR5634762", "SRR5630300", "SRR5630114", "SRR5629890", "SRR5634699", "SRR5630006", "SRR5634715", "SRR5629930", "SRR5630185", "SRR5629904",
                                "SRR5630228", "SRR5630136", "SRR5634385", "SRR5634354", "SRR5630038", "SRR5630237", "SRR5630043", "SRR5634471", "SRR5630382", "SRR5630223", "SRR5629927", "SRR5629819", "SRR5634584", "SRR5634721", "SRR5634693", "SRR5630284", "SRR5634376", "SRR5634566", "SRR5629975", "SRR5630387", "SRR5629854", "SRR5630289", "SRR5629993", "SRR5629961", "SRR5634654", "SRR5630072", "SRR5630172", "SRR5630452", "SRR5630378", "SRR5629785", "SRR5630156", "SRR5634680", "SRR5630390", "SRR5630165", "SRR5630434", "SRR5634606", "SRR5630360", "SRR5629945", "SRR5634695", "SRR5630197", "SRR5630176", "SRR5630407",
                                "SRR5634724", "SRR5630427", "SRR5630221", "SRR5630014", "SRR5630312", "SRR5630069", "SRR5629799", "SRR5630116", "SRR5630142", "SRR5634651", "SRR5634625", "SRR5630417", "SRR5629846", "SRR5634768", "SRR5629917", "SRR5629949", "SRR5634433", "SRR5630216", "SRR5634737", "SRR5634683", "SRR5630117", "SRR5634610", "SRR5629781", "SRR5630372", "SRR5634571", "SRR5634653", "SRR5630240", "SRR5629914", "SRR5634451", "SRR5629888", "SRR5630386", "SRR5630166", "SRR5630297", "SRR5630222", "SRR5630262", "SRR5629937", "SRR5634462", "SRR5630099", "SRR5630206", "SRR5629940", "SRR5630054", "SRR5630027",
                                "SRR5630435", "SRR5634514", "SRR5634722", "SRR5634529", "SRR5634555", "SRR5630380", "SRR5630239", "SRR5629943", "SRR5629966", "SRR5630143", "SRR5634475", "SRR5629946", "SRR5629919", "SRR5634356", "SRR5630229", "SRR5634583", "SRR5634427", "SRR5630182", "SRR5634614", "SRR5629880", "SRR5630283", "SRR5630252", "SRR5630188", "SRR5634805", "SRR5634609", "SRR5630311", "SRR5630269", "SRR5629848", "SRR5630051", "SRR5634756", "SRR5630341", "SRR5634800", "SRR5629906", "SRR5634450", "SRR5630388", "SRR5634537", "SRR5634550", "SRR5634561", "SRR5629855", "SRR5630449", "SRR5629885", "SRR5634630",
                                "SRR5629916", "SRR5629821", "SRR5630144", "SRR5630323", "SRR5634552", "SRR5629988", "SRR5630416", "SRR5634435", "SRR5630272", "SRR5629958", "SRR5630287", "SRR5629922", "SRR5630245", "SRR5634807", "SRR5630044", "SRR5634413", "SRR5630135", "SRR5630203", "SRR5629836", "SRR5634803", "SRR5634748", "SRR5630003", "SRR5634727", "SRR5630247", "SRR5634817", "SRR5630338", "SRR5634374", "SRR5629803", "SRR5630426", "SRR5634635", "SRR5630327", "SRR5634509", "SRR5630062", "SRR5629881", "SRR5634622", "SRR5630091", "SRR5634524", "SRR5630139", "SRR5634812", "SRR5634517", "SRR5630095", "SRR5629925",
                                "SRR5630431", "SRR5634760", "SRR5634434", "SRR5630421", "SRR5630354", "SRR5634637", "SRR5634675", "SRR5630428", "SRR5630004", "SRR5629814", "SRR5630371", "SRR5634506", "SRR5634597", "SRR5629843", "SRR5634458", "SRR5634820", "SRR5634627", "SRR5629923", "SRR5630029", "SRR5630325", "SRR5630405", "SRR5634397", "SRR5629877", "SRR5630149", "SRR5634620", "SRR5630344", "SRR5629983", "SRR5629844", "SRR5634501", "SRR5630021", "SRR5630000", "SRR5630281", "SRR5634811", "SRR5629828", "SRR5630084", "SRR5634608", "SRR5634392", "SRR5629991", "SRR5634763", "SRR5634711", "SRR5634393", "SRR5630432",
                                "SRR5630212", "SRR5634559", "SRR5630308", "SRR5634703", "SRR5634582", "SRR5630010", "SRR5634691", "SRR5629850", "SRR5629956", "SRR5634679", "SRR5630126", "SRR5634769", "SRR5629942", "SRR5634782", "SRR5630196", "SRR5630347", "SRR5629842", "SRR5630352", "SRR5630210", "SRR5630145", "SRR5630296", "SRR5630022", "SRR5630363", "SRR5634496", "SRR5630302", "SRR5634569", "SRR5629812", "SRR5634429", "SRR5630195", "SRR5630007", "SRR5630097", "SRR5630181", "SRR5630412", "SRR5630220", "SRR5629995", "SRR5629899", "SRR5630201", "SRR5629817", "SRR5634432", "SRR5630430", "SRR5634556", "SRR5634771",
                                "SRR5630398", "SRR5629931", "SRR5630446", "SRR5634663", "SRR5634821", "SRR5629809", "SRR5634758", "SRR5630282", "SRR5634565", "SRR5634528", "SRR5629960", "SRR5634741", "SRR5630373", "SRR5630215", "SRR5630305", "SRR5630249", "SRR5630158", "SRR5630274", "SRR5629790", "SRR5634573", "SRR5634564", "SRR5629856", "SRR5629869", "SRR5629816", "SRR5634690", "SRR5629907", "SRR5630049", "SRR5629909", "SRR5629928", "SRR5634358", "SRR5634523", "SRR5634636", "SRR5630019", "SRR5634595", "SRR5629954", "SRR5629795", "SRR5629813", "SRR5634420", "SRR5629957", "SRR5634770", "SRR5630015", "SRR5630108",
                                "SRR5630343", "SRR5634799", "SRR5630140", "SRR5630307", "SRR5629887", "SRR5634361", "SRR5634440", "SRR5630107", "SRR5630092", "SRR5630392", "SRR5629910", "SRR5634407", "SRR5634531", "SRR5630420", "SRR5634692", "SRR5634641", "SRR5630024", "SRR5629793", "SRR5634672", "SRR5634525", "SRR5630442", "SRR5629918", "SRR5634467", "SRR5630001", "SRR5634488", "SRR5634377", "SRR5634776", "SRR5634689", "SRR5630096", "SRR5630133", "SRR5630400", "SRR5629883", "SRR5630298", "SRR5630110", "SRR5629953", "SRR5630189", "SRR5634685", "SRR5630332", "SRR5630253", "SRR5630213", "SRR5634503", "SRR5630187",
                                "SRR5634479", "SRR5629980", "SRR5630391", "SRR5634480", "SRR5630451", "SRR5634439", "SRR5629804", "SRR5629889", "SRR5630234", "SRR5630217", "SRR5630235", "SRR5634700", "SRR5634794", "SRR5634455", "SRR5634619", "SRR5630339", "SRR5630377", "SRR5629818", "SRR5634568", "SRR5630306", "SRR5630005", "SRR5634718", "SRR5634780", "SRR5634394", "SRR5630277", "SRR5634510", "SRR5634381", "SRR5629853", "SRR5629798", "SRR5634621", "SRR5634382", "SRR5630218", "SRR5630288", "SRR5634604", "SRR5634728", "SRR5630279", "SRR5634359", "SRR5630002", "SRR5630443", "SRR5634538", "SRR5629948", "SRR5630134",
                                "SRR5629811", "SRR5629992", "SRR5634472", "SRR5634447", "SRR5634682", "SRR5634370", "SRR5634647", "SRR5634701", "SRR5634493", "SRR5634468", "SRR5634605", "SRR5634572", "SRR5630098", "SRR5634749", "SRR5630319", "SRR5629939", "SRR5630316", "SRR5634624", "SRR5630367"))

# prepare heatmap
p1 = ggplot(cov_file) +
  geom_tile(data = subset(cov_file, cov_file$V7 == "A"),
            aes(x = V3, y = V1, fill = V6)) +
  scale_fill_gradientn(na.value = "grey85", limits = c(0, 2),
                       colours = c("#053061", "#0F437B", "#195696", "#2369AD", "#2F79B5", "#3B89BE", "#4E9AC6", "#6AACD0", "#86BDDA", "#9FCBE1", "#B6D7E8", "#CCE2EE", "#DBEAF2", "#E9F0F4", "#F7F7F7", "#F9ECE5", "#FBE3D4", "#FCD7C2", "#F9C3A9", "#F5B090", "#EF9B7A", "#E58267", "#DA6954", "#CE5045", "#C13639", "#B41D2D", "#9C1127", "#810823", "#67001F"),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Maltose growth",
       x = " ",
       y = "Segregants") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.position = "left",
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

p2 = ggplot(cov_file) +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "R"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "red", limits = c(0, 2),
                       colours = c("#FFFFCC", "#FFF8BC", "#FFF1AC", "#FEEB9C", "#FEE38C", "#FEDC7D", "#FED16E", "#FEC35F", "#FEB54F", "#FDA747", "#FD9A41", "#FD8D3C", "#FC7635", "#FC5F2E", "#F94928", "#F03623", "#E7231E", "#DC151D", "#CE0B21", "#C00225", "#AC0026", "#960026", "#800026"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks = c(0, 1, 2)) +
              new_scale_fill() +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "Y"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "black", limits = c(0, 2),
                       colours = c("white", "grey95", "grey90", "grey85", "grey80", "grey75", "grey70", "grey65", "grey60", "grey55", "grey50", "grey45", "grey40", "grey35", "grey30", "grey25", "grey20", "grey15", "grey10", "black"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Coverage",
       x = " ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

# plot heapmap
grid.arrange(p1, p2, nrow = 1,
             layout_matrix = rbind(c(1, 2)),
             widths = c(0.12, 0.90))

#prepare and plot correlation plot
gen_Y = subset(correlation_file, gen == "Y")
gen_R = subset(correlation_file, gen == "R")
pvalue = t.test(gen_Y$Maltose_growth, gen_R$Maltose_growth)$p.value

ggplot(correlation_file, aes(x = correlation_file$IMA1_cov, y = correlation_file$Maltose_growth)) +
  geom_point(aes(colour = gen), shape = 16, size = 6, alpha = 0.35) +
  geom_errorbar(aes(x = mean(gen_Y$IMA1_cov), ymin = mean(gen_Y$Maltose_growth)-sd(gen_Y$Maltose_growth),
                    ymax = mean(gen_Y$Maltose_growth)+sd(gen_Y$Maltose_growth)),
                width = 0.075, size = 0.75) +
  geom_point(aes(x = mean(gen_Y$IMA1_cov), y = mean(gen_Y$Maltose_growth)),
             shape = 21, size = 6, fill = "white") +
  geom_errorbar(aes(x = mean(gen_R$IMA1_cov), ymin = mean(gen_R$Maltose_growth)-sd(gen_R$Maltose_growth),
                    ymax = mean(gen_R$Maltose_growth)+sd(gen_R$Maltose_growth)),
                width = 0.075, size = 0.75) +
  geom_point(aes(x = mean(gen_R$IMA1_cov), y = mean(gen_R$Maltose_growth)),
             shape = 21, size = 6, fill = "white") +
  geom_rect(aes(xmin = mean(gen_Y$IMA1_cov), xmax = mean(gen_R$IMA1_cov),
                ymin = 1.625, ymax = 1.635)) +
  scale_color_manual(values = c("#FDA747","grey50")) +
  coord_cartesian(xlim = c(-0.1, 1.6), ylim = c(0, 2.1), expand = FALSE) +
  labs(x = "Norm. coverage",
       y = "Maltose growth",
       colour = "Genotype") +
  theme(axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(10, "points"),
        axis.title = element_text(size = 48),
        axis.text = element_text(size = 36),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(size = 1.25, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank()) +
  annotate("Text", x = 0.5, y = 1.65, label = "***", size = 28)

# IMA1 histogram
MALgenes$IMA1_ALLELE = factor(MALgenes$IMA1_ALLELE, levels=c("R", "Y"))
mu = ddply(MALgenes, "IMA1_ALLELE", summarise,
           grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x = Maltose_growth, fill = IMA1_ALLELE, color =
                       IMA1_ALLELE)) +
  geom_histogram(bins = 20, position = "identity", alpha = 0.5) +
  geom_vline(data = mu, aes(xintercept = grp.mean), linetype = "dashed",
             colour = c("#FDA747", "grey50"), size = 1.5) +
  scale_fill_manual(values = c("#FDA747", "grey50")) +
  scale_color_manual(values = c("#FDA747", "grey50")) +
  guides(colour = FALSE) +
  ylim(0, 127) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Maltose growth",
       y = "Number of segregants",
       fill = "Genotype") +
  theme(axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(10, "points"),
        axis.title = element_text(size = 48),
        axis.text = element_text(size = 36),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(size = 1.25, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank())
```



![](./Pics_Maltose/S288C.IMA1cov.heatmap.YR.png)









### RM11 Reference

On the previous analysis we got hints that the IMA1 locus may have been duplicated on RM11. We want to have a closer look at that and we will repeat the analysis using instead RM11 genome as reference.

Here we select ~23 kb: CH408044:1,052,000-1,074,767 (in the EMBL genome would be CH408044:1-22800, somehow me and Ping are working on the reverse complement of the same reference)



```bash
# reads aligmment
while read line; do
	bwa mem -t 72 -K 100000000 \
		../00_REF_genomes/RM11.fa \
		../00_reads/"${line}".R1.tr.fq.gz ../00_reads/"${line}".R2.tr.fq.gz > "${line}".align.sam;
		samtools view -@ 72 -Sb "${line}".align.sam > "${line}".align.bam;
		samtools sort -@ 72 "${line}".align.bam "${line}".align.sort;
	rm *align.sam *align.bam;
done < ../31_IMA1_cov_heatmap/SRR_Acc_List.txt 

bedtools makewindows -g RM11.bedchr -w 100 > RM11.100bp_win.bed
grep CH408044 RM11.100bp_win.bed > RM11.CH408044.100bp_win.bed

for file in *.align.sort.bam; do
	NAME=$(basename $file .align.sort.bam);
	bedtools coverage -a ~/Ping/00_REF_genomes/RM11.CH408044.100bp_win.bed -b $file -mean > $NAME.RM11.100bp_win.cov.bed;
done

grep CH408044 RM11.100bp_win.bed | head -n 228 > RM11.IMA1.locus.100bp_win.bed

# calculate chr7 coverage and normalize windows coverages
for file in *.100bp_win.cov.bed; do
	NAME=$(basename $file .RM11.100bp_win.cov.bed);
	AVG=$(awk '{ total += $4; count++ } END { print total/count }' $file);
	while read line ; do
		VALUE=$(echo $line | cut -f 4 -d ' ');
		NORM=$( echo "scale=8; $VALUE/$AVG" | bc);
		echo $line | tr ' ' '\t' | sed "s/^/$NAME\t/g" | sed "s/$/\t$NORM\tB/g";
	done < <(bedtools intersect -wb -a RM11.IMA1.win.bed -b $file | cut -f 4-8);
done > temp1

while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d ' ');
	GEN=$(echo $line | cut -f 2 -d ' ');
	grep "${SAMPLE}" temp1 | sed "s/$/\t$GEN/g";
done < sample2genotype.tab > temp2

# split dataset and add alpha values
while read line; do
	FIRST=$(echo $line | cut -f 1-4 -d ' ' | tr ' ' '\t');
	COV=$(echo $line | cut -f 5 -d ' ');
	NORM=$(echo $line | cut -f 6 -d ' ');
	GEN=$(echo $line | cut -f 8 -d ' ');
	if [ $GEN == "Y" ]; then
		echo $line$'\t'1 | tr ' ' '\t';
		echo $FIRST$'\t'0.0000000$'\t'0$'\t'B$'\t'R$'\t'0 | tr ' ' '\t';
	elif [ $GEN == "R" ]; then
		echo $line$'\t'1 | tr ' ' '\t';
		echo $FIRST$'\t'0.0000000$'\t'0$'\t'B$'\t'Y$'\t'0 | tr ' ' '\t';
	fi;
done < temp2 | sort -k 8 > temp3

### 2 seconds in python
#	LINE_DB = []
#	with open(args.INFILE) as infile:
#		for line in infile:
#			line = line.strip('\n')
#			STRAIN, CHR, START, STOP, COV, COVAVG, B, GEN = line.split('\t')
#			LINE_DB.append([STRAIN, CHR, str(START), str(STOP), str(COV), str(COVAVG), B, GEN])
#	
#	# open coverage files
#	for ENTRY in LINE_DB:
#		if ENTRY[7] == "R":
#			print('\t'.join(['\t'.join(ENTRY), str(1)]))
#			print('\t'.join([ENTRY[0], ENTRY[1], ENTRY[2], ENTRY[3], "0.0000000", str(0), "B", "Y", str(0)]))
#		if ENTRY[7] == "Y":
#			print('\t'.join(['\t'.join(ENTRY), str(1)]))
#			print('\t'.join([ENTRY[0], ENTRY[1], ENTRY[2], ENTRY[3], "0.0000000", str(0), "B", "R", str(0)]))


cat temp3 ../31_IMA1_cov_heatmap/02_coverage/allsamples.100bp.cov.malgrowth.YR.tab > allsamples.100bp.cov.heatmap.YR.tab
rm temp1 temp2 temp3

# prepare correlation plot
echo CH408044$'\t'9100$'\t'20500 > RM11.IMA1.shortallele.bed
for file in 02_coverage/*.100bp_win.cov.bed; do
	NAME=$(basename $file .RM11.100bp_win.cov.bed);
	AVG=$(awk '{ total += $4; count++ } END { print total/count }' $file);
	while read line; do
		VALUE=$(echo $line | cut -f 4 -d ' ');
		NORM=$( echo "scale=8; $VALUE/$AVG" | bc);
		echo $NORM;
	done < <(bedtools intersect -wb -a RM11.IMA1.shortallele.bed -b $file | cut -f 4-8) |\
    	awk '{ total += $1; count++ } END { print total/count }' - |\
        sed "s/^/$NAME\t/g";
done > allsamples.RM11.IMA1.shortallele.YR.tab.1

while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d " ");
	GEN=$(echo $line | cut -f 2 -d " ");
	grep "${SAMPLE}" allsamples.RM11.IMA1.shortallele.YR.tab.1 | sed "s/$/\t$GEN/g";
done < sample2genotype.tab > allsamples.RM11.IMA1.shortallele.YR.tab

echo Maltose_growth$'\t'IMA1_cov$'\t'gen > allsamples.100bp.cov.correlation.YR.tab;
while read line; do
	SAMPLE=$(echo $line | cut -f 1 -d " ");
	PHEN=$(echo $line | cut -f 2 -d " ");
	grep "${SAMPLE}" allsamples.RM11.IMA1.shortallele.YR.tab | sed "s/$SAMPLE/$PHEN/g";
done < ../31_IMA1_cov_heatmap/02_coveragesample2malgrowth.tab >> allsamples.100bp.cov.correlation.YR.tab
```





**RM11.IMA1cov.heatmap.R**

```R
library("dplyr")
library("ggplot2")
library("ggpubr")
library("grid")
library("gridExtra")
library("plyr")
library("scales")
library("RColorBrewer")
library("ggnewscale")
library("tidyverse")

# upload and organize
setwd("/media/DISK2-3TB//Ping/32_IMA1_cov_heatmap_RM11//")
cov_file = read.delim("allsamples.100bp.cov.heatmap.YR.tab", header = FALSE)
correlation_file = read.delim("allsamples.100bp.cov.correlation.YR.tab", header = TRUE)
MALgenes = read.delim("strains_MAL_gen2phe.txt", header = TRUE)

cov_file$V1 = factor(cov_file$V1,
                     levels = c("SRR5634755", "SRR5634774", "SRR5630056", "SRR5634576", "SRR5630345", "SRR5634504", "SRR5630053", "SRR5630190", "SRR5630087", "SRR5634745", "SRR5629932", "SRR5634542", "SRR5630013", "SRR5634491", "SRR5634513", "SRR5634546", "SRR5634530", "SRR5630171", "SRR5634612", "SRR5630175", "SRR5634823", "SRR5630349", "SRR5634650", "SRR5634752", "SRR5629824", "SRR5634416", "SRR5634482", "SRR5629863", "SRR5629972", "SRR5634730", "SRR5630447", "SRR5634498", "SRR5630208", "SRR5634441", "SRR5630164", "SRR5634349", "SRR5634368", "SRR5634725", "SRR5634767", "SRR5629868", "SRR5634391", "SRR5634638",
                                "SRR5629859", "SRR5630263", "SRR5630419", "SRR5634667", "SRR5634507", "SRR5634666", "SRR5630353", "SRR5634460", "SRR5629879", "SRR5629872", "SRR5630424", "SRR5634547", "SRR5634657", "SRR5630200", "SRR5634759", "SRR5634661", "SRR5630362", "SRR5634676", "SRR5629810", "SRR5630147", "SRR5630299", "SRR5629845", "SRR5630331", "SRR5630356", "SRR5630399", "SRR5634466", "SRR5630052", "SRR5630103", "SRR5630291", "SRR5634408", "SRR5630093", "SRR5629950", "SRR5634532", "SRR5629838", "SRR5630368", "SRR5630076", "SRR5634463", "SRR5630211", "SRR5630040", "SRR5629851", "SRR5634421", "SRR5634645",
                                "SRR5634686", "SRR5629801", "SRR5634586", "SRR5630121", "SRR5630034", "SRR5629835", "SRR5630351", "SRR5630023", "SRR5629830", "SRR5630032", "SRR5630444", "SRR5634431", "SRR5630438", "SRR5634596", "SRR5630320", "SRR5634364", "SRR5634459", "SRR5634540", "SRR5630048", "SRR5634786", "SRR5634558", "SRR5629805", "SRR5630192", "SRR5630074", "SRR5630148", "SRR5630105", "SRR5629891", "SRR5634362", "SRR5629920", "SRR5634438", "SRR5634406", "SRR5629864", "SRR5634422", "SRR5630425", "SRR5630028", "SRR5630326", "SRR5634499", "SRR5634706", "SRR5630170", "SRR5630179", "SRR5630233", "SRR5634810",
                                "SRR5630128", "SRR5630225", "SRR5629794", "SRR5630243", "SRR5630155", "SRR5630384", "SRR5630264", "SRR5634570", "SRR5630045", "SRR5629959", "SRR5630374", "SRR5630395", "SRR5634743", "SRR5630246", "SRR5629934", "SRR5634601", "SRR5634402", "SRR5629973", "SRR5634414", "SRR5630129", "SRR5634673", "SRR5634519", "SRR5630276", "SRR5629981", "SRR5634696", "SRR5634740", "SRR5630254", "SRR5630337", "SRR5634369", "SRR5634618", "SRR5629829", "SRR5630304", "SRR5634575", "SRR5634735", "SRR5634640", "SRR5629987", "SRR5634516", "SRR5634736", "SRR5629825", "SRR5629978", "SRR5630184", "SRR5634409",
                                "SRR5634792", "SRR5629791", "SRR5629871", "SRR5634581", "SRR5634616", "SRR5634452", "SRR5634723", "SRR5634512", "SRR5634642", "SRR5630266", "SRR5630448", "SRR5634588", "SRR5634430", "SRR5629874", "SRR5634813", "SRR5634500", "SRR5630280", "SRR5630379", "SRR5634753", "SRR5630270", "SRR5629849", "SRR5630008", "SRR5630375", "SRR5629913", "SRR5634824", "SRR5634623", "SRR5630086", "SRR5630396", "SRR5630365", "SRR5634388", "SRR5630075", "SRR5634567", "SRR5634802", "SRR5634423", "SRR5629941", "SRR5634660", "SRR5634350", "SRR5634474", "SRR5630146", "SRR5630389", "SRR5630415", "SRR5634481",
                                "SRR5634634", "SRR5630242", "SRR5630406", "SRR5629837", "SRR5630257", "SRR5630031", "SRR5630065", "SRR5630322", "SRR5634505", "SRR5630080", "SRR5634539", "SRR5634456", "SRR5630376", "SRR5634415", "SRR5634665", "SRR5630113", "SRR5634476", "SRR5630067", "SRR5629806", "SRR5629886", "SRR5630361", "SRR5630039", "SRR5634494", "SRR5629847", "SRR5629788", "SRR5634424", "SRR5630324", "SRR5630369", "SRR5630259", "SRR5630329", "SRR5634502", "SRR5630012", "SRR5634508", "SRR5634681", "SRR5634777", "SRR5634744", "SRR5630152", "SRR5634687", "SRR5630081", "SRR5634708", "SRR5634535", "SRR5634428",
                                "SRR5634348", "SRR5634716", "SRR5630255", "SRR5630385", "SRR5634822", "SRR5629903", "SRR5634464", "SRR5634403", "SRR5630068", "SRR5634399", "SRR5630109", "SRR5630408", "SRR5634677", "SRR5629968", "SRR5634798", "SRR5630106", "SRR5630073", "SRR5634515", "SRR5630290", "SRR5629971", "SRR5629974", "SRR5629858", "SRR5629979", "SRR5630104", "SRR5634603", "SRR5634492", "SRR5634707", "SRR5630413", "SRR5634384", "SRR5634746", "SRR5629841", "SRR5634779", "SRR5629797", "SRR5629834", "SRR5634347", "SRR5629820", "SRR5630046", "SRR5630295", "SRR5629840", "SRR5630241", "SRR5630205", "SRR5630265",
                                "SRR5634655", "SRR5630357", "SRR5630077", "SRR5634375", "SRR5629815", "SRR5630178", "SRR5629792", "SRR5630101", "SRR5630433", "SRR5630401", "SRR5630403", "SRR5630286", "SRR5630303", "SRR5634719", "SRR5634549", "SRR5634790", "SRR5634554", "SRR5629969", "SRR5634527", "SRR5630102", "SRR5630364", "SRR5634437", "SRR5634580", "SRR5630199", "SRR5634478", "SRR5634643", "SRR5634592", "SRR5634607", "SRR5630273", "SRR5629857", "SRR5630227", "SRR5629782", "SRR5630333", "SRR5629921", "SRR5630334", "SRR5629867", "SRR5630082", "SRR5630232", "SRR5630397", "SRR5629962", "SRR5634713", "SRR5634383",
                                "SRR5630250", "SRR5629882", "SRR5630318", "SRR5634602", "SRR5634633", "SRR5630026", "SRR5630226", "SRR5630132", "SRR5630383", "SRR5634649", "SRR5629823", "SRR5634521", "SRR5634386", "SRR5630177", "SRR5629802", "SRR5634714", "SRR5629951", "SRR5634405", "SRR5629915", "SRR5634511", "SRR5634652", "SRR5634754", "SRR5634747", "SRR5630301", "SRR5630244", "SRR5629964", "SRR5634738", "SRR5630358", "SRR5634390", "SRR5630180", "SRR5629832", "SRR5634750", "SRR5634709", "SRR5630214", "SRR5630278", "SRR5634670", "SRR5629929", "SRR5630268", "SRR5630439", "SRR5634731", "SRR5629911", "SRR5634662",
                                "SRR5634401", "SRR5634587", "SRR5630251", "SRR5630445", "SRR5629984", "SRR5634366", "SRR5630112", "SRR5630138", "SRR5634598", "SRR5630411", "SRR5634632", "SRR5634473", "SRR5634704", "SRR5630123", "SRR5634694", "SRR5630169", "SRR5630154", "SRR5634804", "SRR5630159", "SRR5634659", "SRR5629908", "SRR5630118", "SRR5630017", "SRR5634577", "SRR5629963", "SRR5629822", "SRR5629997", "SRR5634639", "SRR5630441", "SRR5630094", "SRR5630410", "SRR5634664", "SRR5630248", "SRR5630370", "SRR5634674", "SRR5634668", "SRR5630328", "SRR5634720", "SRR5634726", "SRR5630224", "SRR5629831", "SRR5634520",
                                "SRR5629783", "SRR5634626", "SRR5630186", "SRR5630089", "SRR5634454", "SRR5634648", "SRR5634465", "SRR5634372", "SRR5630409", "SRR5630137", "SRR5634757", "SRR5634772", "SRR5630041", "SRR5630124", "SRR5634352", "SRR5629875", "SRR5634742", "SRR5634784", "SRR5630059", "SRR5630209", "SRR5634449", "SRR5634410", "SRR5630422", "SRR5630258", "SRR5634543", "SRR5634773", "SRR5629839", "SRR5634453", "SRR5630236", "SRR5634443", "SRR5629796", "SRR5630191", "SRR5634710", "SRR5629976", "SRR5634551", "SRR5630204", "SRR5630058", "SRR5630267", "SRR5629860", "SRR5634712", "SRR5630157", "SRR5634404",
                                "SRR5630231", "SRR5634815", "SRR5630294", "SRR5634379", "SRR5634518", "SRR5634669", "SRR5634761", "SRR5630055", "SRR5630256", "SRR5629938", "SRR5634729", "SRR5634688", "SRR5630011", "SRR5629826", "SRR5634495", "SRR5630437", "SRR5630153", "SRR5630033", "SRR5634658", "SRR5634411", "SRR5634387", "SRR5634426", "SRR5630070", "SRR5630293", "SRR5634796", "SRR5630330", "SRR5630079", "SRR5629989", "SRR5630402", "SRR5629878", "SRR5629936", "SRR5630090", "SRR5634764", "SRR5630078", "SRR5634357", "SRR5630125", "SRR5630335", "SRR5634485", "SRR5630063", "SRR5634628", "SRR5630271", "SRR5634826",
                                "SRR5629967", "SRR5634819", "SRR5629862", "SRR5630141", "SRR5634733", "SRR5629970", "SRR5629827", "SRR5630061", "SRR5629998", "SRR5630366", "SRR5630016", "SRR5630161", "SRR5630150", "SRR5634442", "SRR5630085", "SRR5634395", "SRR5634617", "SRR5630194", "SRR5630450", "SRR5634544", "SRR5634563", "SRR5634351", "SRR5629896", "SRR5629861", "SRR5629999", "SRR5629870", "SRR5634425", "SRR5634717", "SRR5634522", "SRR5630115", "SRR5630359", "SRR5634553", "SRR5634470", "SRR5630163", "SRR5634806", "SRR5630219", "SRR5630066", "SRR5630037", "SRR5630202", "SRR5634739", "SRR5634446", "SRR5634613",
                                "SRR5630285", "SRR5630350", "SRR5630183", "SRR5634365", "SRR5630160", "SRR5634705", "SRR5630394", "SRR5630310", "SRR5630088", "SRR5629866", "SRR5634417", "SRR5634809", "SRR5634797", "SRR5634615", "SRR5634698", "SRR5629994", "SRR5630313", "SRR5630292", "SRR5634487", "SRR5634367", "SRR5629892", "SRR5634557", "SRR5630393", "SRR5629808", "SRR5634461", "SRR5630340", "SRR5629865", "SRR5634600", "SRR5630346", "SRR5630317", "SRR5634629", "SRR5629807", "SRR5629965", "SRR5629986", "SRR5630275", "SRR5634778", "SRR5629947", "SRR5629996", "SRR5630314", "SRR5629977", "SRR5629955", "SRR5634702",
                                "SRR5630309", "SRR5630440", "SRR5629894", "SRR5629884", "SRR5634548", "SRR5630381", "SRR5634545", "SRR5634419", "SRR5634793", "SRR5634541", "SRR5630151", "SRR5630064", "SRR5630230", "SRR5634765", "SRR5634671", "SRR5634489", "SRR5634560", "SRR5634477", "SRR5630047", "SRR5634497", "SRR5634534", "SRR5630083", "SRR5630127", "SRR5634801", "SRR5634398", "SRR5634825", "SRR5630174", "SRR5634469", "SRR5634816", "SRR5634644", "SRR5634611", "SRR5630238", "SRR5634814", "SRR5630261", "SRR5634631", "SRR5629990", "SRR5629924", "SRR5630418", "SRR5629786", "SRR5630050", "SRR5629933", "SRR5629926",
                                "SRR5630036", "SRR5629833", "SRR5634355", "SRR5634536", "SRR5634400", "SRR5634766", "SRR5629985", "SRR5630018", "SRR5630100", "SRR5634396", "SRR5630193", "SRR5634363", "SRR5629873", "SRR5634599", "SRR5634562", "SRR5634808", "SRR5634697", "SRR5634732", "SRR5634389", "SRR5634734", "SRR5630071", "SRR5630404", "SRR5630042", "SRR5630414", "SRR5630057", "SRR5634353", "SRR5630020", "SRR5634589", "SRR5630315", "SRR5629905", "SRR5634486", "SRR5634371", "SRR5634762", "SRR5630300", "SRR5630114", "SRR5629890", "SRR5634699", "SRR5630006", "SRR5634715", "SRR5629930", "SRR5630185", "SRR5629904",
                                "SRR5630228", "SRR5630136", "SRR5634385", "SRR5634354", "SRR5630038", "SRR5630237", "SRR5630043", "SRR5634471", "SRR5630382", "SRR5630223", "SRR5629927", "SRR5629819", "SRR5634584", "SRR5634721", "SRR5634693", "SRR5630284", "SRR5634376", "SRR5634566", "SRR5629975", "SRR5630387", "SRR5629854", "SRR5630289", "SRR5629993", "SRR5629961", "SRR5634654", "SRR5630072", "SRR5630172", "SRR5630452", "SRR5630378", "SRR5629785", "SRR5630156", "SRR5634680", "SRR5630390", "SRR5630165", "SRR5630434", "SRR5634606", "SRR5630360", "SRR5629945", "SRR5634695", "SRR5630197", "SRR5630176", "SRR5630407",
                                "SRR5634724", "SRR5630427", "SRR5630221", "SRR5630014", "SRR5630312", "SRR5630069", "SRR5629799", "SRR5630116", "SRR5630142", "SRR5634651", "SRR5634625", "SRR5630417", "SRR5629846", "SRR5634768", "SRR5629917", "SRR5629949", "SRR5634433", "SRR5630216", "SRR5634737", "SRR5634683", "SRR5630117", "SRR5634610", "SRR5629781", "SRR5630372", "SRR5634571", "SRR5634653", "SRR5630240", "SRR5629914", "SRR5634451", "SRR5629888", "SRR5630386", "SRR5630166", "SRR5630297", "SRR5630222", "SRR5630262", "SRR5629937", "SRR5634462", "SRR5630099", "SRR5630206", "SRR5629940", "SRR5630054", "SRR5630027",
                                "SRR5630435", "SRR5634514", "SRR5634722", "SRR5634529", "SRR5634555", "SRR5630380", "SRR5630239", "SRR5629943", "SRR5629966", "SRR5630143", "SRR5634475", "SRR5629946", "SRR5629919", "SRR5634356", "SRR5630229", "SRR5634583", "SRR5634427", "SRR5630182", "SRR5634614", "SRR5629880", "SRR5630283", "SRR5630252", "SRR5630188", "SRR5634805", "SRR5634609", "SRR5630311", "SRR5630269", "SRR5629848", "SRR5630051", "SRR5634756", "SRR5630341", "SRR5634800", "SRR5629906", "SRR5634450", "SRR5630388", "SRR5634537", "SRR5634550", "SRR5634561", "SRR5629855", "SRR5630449", "SRR5629885", "SRR5634630",
                                "SRR5629916", "SRR5629821", "SRR5630144", "SRR5630323", "SRR5634552", "SRR5629988", "SRR5630416", "SRR5634435", "SRR5630272", "SRR5629958", "SRR5630287", "SRR5629922", "SRR5630245", "SRR5634807", "SRR5630044", "SRR5634413", "SRR5630135", "SRR5630203", "SRR5629836", "SRR5634803", "SRR5634748", "SRR5630003", "SRR5634727", "SRR5630247", "SRR5634817", "SRR5630338", "SRR5634374", "SRR5629803", "SRR5630426", "SRR5634635", "SRR5630327", "SRR5634509", "SRR5630062", "SRR5629881", "SRR5634622", "SRR5630091", "SRR5634524", "SRR5630139", "SRR5634812", "SRR5634517", "SRR5630095", "SRR5629925",
                                "SRR5630431", "SRR5634760", "SRR5634434", "SRR5630421", "SRR5630354", "SRR5634637", "SRR5634675", "SRR5630428", "SRR5630004", "SRR5629814", "SRR5630371", "SRR5634506", "SRR5634597", "SRR5629843", "SRR5634458", "SRR5634820", "SRR5634627", "SRR5629923", "SRR5630029", "SRR5630325", "SRR5630405", "SRR5634397", "SRR5629877", "SRR5630149", "SRR5634620", "SRR5630344", "SRR5629983", "SRR5629844", "SRR5634501", "SRR5630021", "SRR5630000", "SRR5630281", "SRR5634811", "SRR5629828", "SRR5630084", "SRR5634608", "SRR5634392", "SRR5629991", "SRR5634763", "SRR5634711", "SRR5634393", "SRR5630432",
                                "SRR5630212", "SRR5634559", "SRR5630308", "SRR5634703", "SRR5634582", "SRR5630010", "SRR5634691", "SRR5629850", "SRR5629956", "SRR5634679", "SRR5630126", "SRR5634769", "SRR5629942", "SRR5634782", "SRR5630196", "SRR5630347", "SRR5629842", "SRR5630352", "SRR5630210", "SRR5630145", "SRR5630296", "SRR5630022", "SRR5630363", "SRR5634496", "SRR5630302", "SRR5634569", "SRR5629812", "SRR5634429", "SRR5630195", "SRR5630007", "SRR5630097", "SRR5630181", "SRR5630412", "SRR5630220", "SRR5629995", "SRR5629899", "SRR5630201", "SRR5629817", "SRR5634432", "SRR5630430", "SRR5634556", "SRR5634771",
                                "SRR5630398", "SRR5629931", "SRR5630446", "SRR5634663", "SRR5634821", "SRR5629809", "SRR5634758", "SRR5630282", "SRR5634565", "SRR5634528", "SRR5629960", "SRR5634741", "SRR5630373", "SRR5630215", "SRR5630305", "SRR5630249", "SRR5630158", "SRR5630274", "SRR5629790", "SRR5634573", "SRR5634564", "SRR5629856", "SRR5629869", "SRR5629816", "SRR5634690", "SRR5629907", "SRR5630049", "SRR5629909", "SRR5629928", "SRR5634358", "SRR5634523", "SRR5634636", "SRR5630019", "SRR5634595", "SRR5629954", "SRR5629795", "SRR5629813", "SRR5634420", "SRR5629957", "SRR5634770", "SRR5630015", "SRR5630108",
                                "SRR5630343", "SRR5634799", "SRR5630140", "SRR5630307", "SRR5629887", "SRR5634361", "SRR5634440", "SRR5630107", "SRR5630092", "SRR5630392", "SRR5629910", "SRR5634407", "SRR5634531", "SRR5630420", "SRR5634692", "SRR5634641", "SRR5630024", "SRR5629793", "SRR5634672", "SRR5634525", "SRR5630442", "SRR5629918", "SRR5634467", "SRR5630001", "SRR5634488", "SRR5634377", "SRR5634776", "SRR5634689", "SRR5630096", "SRR5630133", "SRR5630400", "SRR5629883", "SRR5630298", "SRR5630110", "SRR5629953", "SRR5630189", "SRR5634685", "SRR5630332", "SRR5630253", "SRR5630213", "SRR5634503", "SRR5630187",
                                "SRR5634479", "SRR5629980", "SRR5630391", "SRR5634480", "SRR5630451", "SRR5634439", "SRR5629804", "SRR5629889", "SRR5630234", "SRR5630217", "SRR5630235", "SRR5634700", "SRR5634794", "SRR5634455", "SRR5634619", "SRR5630339", "SRR5630377", "SRR5629818", "SRR5634568", "SRR5630306", "SRR5630005", "SRR5634718", "SRR5634780", "SRR5634394", "SRR5630277", "SRR5634510", "SRR5634381", "SRR5629853", "SRR5629798", "SRR5634621", "SRR5634382", "SRR5630218", "SRR5630288", "SRR5634604", "SRR5634728", "SRR5630279", "SRR5634359", "SRR5630002", "SRR5630443", "SRR5634538", "SRR5629948", "SRR5630134",
                                "SRR5629811", "SRR5629992", "SRR5634472", "SRR5634447", "SRR5634682", "SRR5634370", "SRR5634647", "SRR5634701", "SRR5634493", "SRR5634468", "SRR5634605", "SRR5634572", "SRR5630098", "SRR5634749", "SRR5630319", "SRR5629939", "SRR5630316", "SRR5634624", "SRR5630367"))

# prepare heatmap
p1 = ggplot(cov_file) +
  geom_tile(data = subset(cov_file, cov_file$V7 == "A"),
            aes(x = V3, y = V1, fill = V6)) +
  scale_fill_gradientn(na.value = "grey85", limits = c(0, 2),
                       colours = c("#053061", "#0F437B", "#195696", "#2369AD", "#2F79B5", "#3B89BE", "#4E9AC6", "#6AACD0", "#86BDDA", "#9FCBE1", "#B6D7E8", "#CCE2EE", "#DBEAF2", "#E9F0F4", "#F7F7F7", "#F9ECE5", "#FBE3D4", "#FCD7C2", "#F9C3A9", "#F5B090", "#EF9B7A", "#E58267", "#DA6954", "#CE5045", "#C13639", "#B41D2D", "#9C1127", "#810823", "#67001F"),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Maltose growth",
       x = " ",
       y = "Segregants") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.position = "left",
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

p2 = ggplot(cov_file) +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "R"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "red", limits = c(0, 2),
                       colours = c("#FFFFCC", "#FFF8BC", "#FFF1AC", "#FEEB9C", "#FEE38C", "#FEDC7D", "#FED16E", "#FEC35F", "#FEB54F", "#FDA747", "#FD9A41", "#FD8D3C", "#FC7635", "#FC5F2E", "#F94928", "#F03623", "#E7231E", "#DC151D", "#CE0B21", "#C00225", "#AC0026", "#960026", "#800026"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks = c(0, 1, 2)) +
              new_scale_fill() +
  geom_tile(data =  subset(cov_file, cov_file$V7 == "B" & cov_file$V8 == "Y"),
            aes(x = V3, y = V1, fill = V6, alpha = V9 )) +
  scale_fill_gradientn(na.value = "black", limits = c(0, 2),
                       colours = c("white", "grey95", "grey90", "grey85", "grey80", "grey75", "grey70", "grey65", "grey60", "grey55", "grey50", "grey45", "grey40", "grey35", "grey30", "grey25", "grey20", "grey15", "grey10", "black"),
                       values = c(0, 0.00000001, 1, 2),
                       breaks=c(0, 1, 2)) +
  scale_x_continuous(labels = comma) +
  scale_x_reverse() +
  coord_cartesian(expand = FALSE) +
  labs(fill = "Coverage",
       x = " ") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA, size = 0.75)) +
  guides(fill = guide_colourbar(ticks.colour = "black",
                                frame.colour = "black",
                                frame.linewidth = 1))

# plot heapmap
grid.arrange(p1, p2, nrow = 1,
             layout_matrix = rbind(c(1, 2)),
             widths = c(0.12, 0.90))

#prepare and plot correlation plot
gen_Y = subset(correlation_file, gen == "Y")
gen_R = subset(correlation_file, gen == "R")
pvalue = t.test(gen_Y$Maltose_growth, gen_R$Maltose_growth)$p.value

ggplot(correlation_file, aes(x = correlation_file$IMA1_cov, y = correlation_file$Maltose_growth)) +
  geom_point(aes(colour = gen), shape = 16, size = 6, alpha = 0.35) +
  geom_errorbar(aes(x = mean(gen_Y$IMA1_cov), ymin = mean(gen_Y$Maltose_growth)-sd(gen_Y$Maltose_growth),
                    ymax = mean(gen_Y$Maltose_growth)+sd(gen_Y$Maltose_growth)),
                width = 0.075, size = 0.75) +
  geom_point(aes(x = mean(gen_Y$IMA1_cov), y = mean(gen_Y$Maltose_growth)),
             shape = 21, size = 6, fill = "white") +
  geom_errorbar(aes(x = mean(gen_R$IMA1_cov), ymin = mean(gen_R$Maltose_growth)-sd(gen_R$Maltose_growth),
                    ymax = mean(gen_R$Maltose_growth)+sd(gen_R$Maltose_growth)),
                width = 0.075, size = 0.75) +
  geom_point(aes(x = mean(gen_R$IMA1_cov), y = mean(gen_R$Maltose_growth)),
             shape = 21, size = 6, fill = "white") +
  geom_rect(aes(xmin = mean(gen_Y$IMA1_cov), xmax = mean(gen_R$IMA1_cov),
                ymin = 1.625, ymax = 1.635)) +
  scale_color_manual(values = c("#FDA747","grey50")) +
  coord_cartesian(xlim = c(-0.1, 1.25), ylim = c(0, 2.1), expand = FALSE) +
  labs(x = "Norm. coverage",
       y = "Maltose growth",
       colour = "Genotype") +
  theme(axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(10, "points"),
        axis.title = element_text(size = 48),
        axis.text = element_text(size = 36),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(size = 1.25, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank()) +
  annotate("Text", x = 0.5, y = 1.65, label = "***", size = 28)

# IMA1 histogram
MALgenes$IMA1_ALLELE = factor(MALgenes$IMA1_ALLELE, levels=c("R", "Y"))
mu = ddply(MALgenes, "IMA1_ALLELE", summarise,
           grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x = Maltose_growth, fill = IMA1_ALLELE, color =
                       IMA1_ALLELE)) +
  geom_histogram(bins = 20, position = "identity", alpha = 0.5) +
  geom_vline(data = mu, aes(xintercept = grp.mean), linetype = "dashed",
             colour = c("#FDA747", "grey50"), size = 1.5) +
  scale_fill_manual(values = c("#FDA747", "grey50")) +
  scale_color_manual(values = c("#FDA747", "grey50")) +
  guides(colour = FALSE) +
  ylim(0, 127) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Maltose growth",
       y = "Number of segregants",
       fill = "Genotype") +
  theme(axis.ticks = element_line(size = 0.75),
        axis.ticks.length = unit(10, "points"),
        axis.title = element_text(size = 48),
        axis.text = element_text(size = 36),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 24),
        legend.direction = "horizontal",
        legend.position = c(0.5, 0.95),
        legend.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(size = 1.25, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_blank())
```





![](./Pics_Maltose/RM11.IMA1cov.heatmap.YR.png)









# _____________________________________________________________________________________________



## CUSTOM SCRIPTS

Hereafter I include all the custom scripts used for all the analyses reported in the same order as they appear in this manual. R scripts to generate the figures stay above the corresponding figure(s), for ease of reproduction. All scripts (R scripts for visualization included) are written by Andrea Del Cortona, unless otherwise specified.



### SelectList_Fasta.pl

```perl
#! /usr/bin/env perl

=head1 Description

  Reads as par1 a proteome (multiple fasta), as par 2 a list (txt file) 
  with entries that will be kept from the dataset.
  
  Result to STDOUT


=cut

#------------------------------------------------------------

sub usage( $ )
  {
    print "$_[0]\n";
    system("pod2text $0");
    exit(1);
  }

#------------------------------------------------------------

(-e "$ARGV[0]" && -e "$ARGV[1]") || &usage("input file\n");
print STDERR "* Reading $ARGV[0]\n";
%prot=&fasta2hash("$ARGV[0]");
print STDERR "* Selecting seqs...\n";
open (LIST, "$ARGV[1]");
@LIST=<LIST>;
close (LIST);
chomp(@LIST);

foreach $keep (@LIST)
 {
# $keep =~ s/[\+-]$//g;
  $keep =~ s/^>//g;
# print STDERR "List $keep\n";
 if (exists $prot{$keep})
  {
  print STDOUT ">$keep\n$prot{$keep}\n";
  $n++;
  }
  else
  {
  #print STDERR "$keep not found\n";
  }
 }
 
print STDERR "* $n saved (".(scalar(@LIST)).") in list.\n";
 
# SUB

#------------------------------------------------------------
# Reads in an entire (multiple) fasta file and returns a hash in which
# the keys are the identifiers of the sequences (without the '>')
# and the values are the sequences themselves. Created by cesim 14/01/2002
#

sub fasta2hash ( $ )
 {
  my ($file,$key,$value);
  my (%fasta_hash);
  $file=$_[0];
  open (IN,$file);
  while (<IN>)
   {
    chomp;
    if (/^>(\S+)/)
     {
      $key=$1;
 #     print STDERR "FH key $key\t"; #debug
     } #if (/^>(\w)$/)
    else
     {
      (defined $key) || die "File $file is not a fasta file!";
      s/\s+//g;
      $fasta_hash{$key}.=uc($_);
     } #else
   } #while (<IN>)
  close IN;
  return (%fasta_hash);
 } #fasta2hash ( $ )
```



### fastq-remove-orphans.pl

```perl
#! /usr/bin/perl
# Victor Amin 2009

use strict;
use warnings;

use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;

my %options;
getopts('1:2:h', \%options);

if ($options{h} || !$options{1} || !$options{2}) {Getopt::Std->version_mess(); HELP_MESSAGE(\*STDERR)}
sub HELP_MESSAGE {
  my $fh = shift;
  print $fh "\nSplit ophaned reads out of a pair of FASTQ files. Counts to STDOUT.\n";
  print $fh "\tOPTIONS:\n";
  print $fh "\t-1 [FASTQ1] [required]\n";
  print $fh "\t-2 [FASTQ2] [required]\n";
  print $fh "\nProperly paired FASTQs are outputted to paired_*, orphans to orphaned_*\n\n";
  exit;
}

open FASTQ1, "<$options{1}" or die "\nThere was a problem opening the FASTQ file: $!\n";
open FASTQ2, "<$options{2}" or die "\nThere was a problem opening the FASTQ file: $!\n";

open PAIRED1, ">paired_$options{1}" or die "\nThere was a problem opening the output file: $!\n";
open PAIRED2, ">paired_$options{2}" or die "\nThere was a problem opening the output file: $!\n";

open ORPHANED1, ">orphaned_$options{1}" or die "\nThere was a problem opening the output file: $!\n";
open ORPHANED2, ">orphaned_$options{2}" or die "\nThere was a problem opening the output file: $!\n";

my $SEQ_MODE = 1;
my $QUAL_MODE = 2;
my $mode = 1;

my $reads_1 = 0;
my $lines = 0;
my $ident;
my $sequence;
my $quality;

my %sequences_1;
my %qualities_1;
print STDERR "\nLoading first FASTQ...\n";
while (<FASTQ1>) {
    chomp;
    if (/^\@/ && $mode == $SEQ_MODE) {
    /\@([^ \/]+?)( |\/)/; # rather than substitute capture AU 16/08/2012
    # chop; this is to remove the 1 or 2 from /1 or /2 if the reads are in that format. These reads are in the format @ident 1:N:0:x or @ident 2:N:0:x AU 16/08/2012
    $ident = $1; # ident = capture AU 16/08/2012
    $reads_1++;
    } elsif (/^\+/) {
    $mode = $QUAL_MODE;
    } elsif ($mode == $SEQ_MODE) {
    $sequence .= $_;
    $lines++;
    } elsif ($mode == $QUAL_MODE) {
    $quality .= $_;
    $lines--;
    if ($lines == 0) {
      $mode = $SEQ_MODE;
      $sequences_1{$ident} = $sequence;
      $qualities_1{$ident} = $quality;
      $sequence = '';
      $quality = '';
    }
    } else {
    die "\nError reading file.\n";
    }
}

my $reads_2 = 0;

my %sequences_2;
my %qualities_2;
print STDERR "\nLoading second FASTQ...\n";
while (<FASTQ2>) {
    chomp;
    if (/^\@/ && $mode == $SEQ_MODE) {
    /\@([^ \/]+?)( |\/)/; # rather than substitute capture AU 16/08/2012
    # chop; this is to remove the 1 or 2 from /1 or /2 if the reads are in that format. These reads are in the format @ident 1:N:0:x or @ident 2:N:0:x AU 16/08/2012
    $ident = $1; # ident = capture AU 16/08/2012
    $reads_2++;
    } elsif (/^\+/) {
    $mode = $QUAL_MODE;
    } elsif ($mode == $SEQ_MODE) {
    $sequence .= $_;
    $lines++;
    } elsif ($mode == $QUAL_MODE) {
    $quality .= $_;
    $lines--;
    if ($lines == 0) {
      $mode = $SEQ_MODE;
      $sequences_2{$ident} = $sequence;
      $qualities_2{$ident} = $quality;
      $sequence = '';
      $quality = '';
    }
    } else {
    die "\nError reading file.\n";
    }
}

my $paired;
print STDERR "\nPrinting paired reads...\n";
for $ident (keys %sequences_1) {
  if (exists $sequences_2{$ident}) {
    print PAIRED1 "\@${ident} 1\n$sequences_1{$ident}\n\+${ident} 1\n$qualities_1{$ident}\n";
    print PAIRED2 "\@${ident} 2\n$sequences_2{$ident}\n\+${ident} 2\n$qualities_2{$ident}\n";
    delete $sequences_1{$ident};
    delete $sequences_2{$ident};
    $paired++;
  }
}

print STDERR "\nPrinting orphaned reads...\n";
my $orphaned_1 = 0;
for $ident (keys %sequences_1) {
  print ORPHANED1  "\@${ident} 1\n$sequences_1{$ident}\n\+${ident} 1\n$qualities_1{$ident}\n";
  $orphaned_1++;
}

my $orphaned_2 = 0;
for $ident (keys %sequences_2) {
  print ORPHANED2  "\@${ident} 2\n$sequences_2{$ident}\n\+${ident} 2\n$qualities_2{$ident}\n";
  $orphaned_2++
}

print "\nReads 1: $reads_1\nOrphans 1: $orphaned_1\nReads 2: $reads_2\nOrphaned 2: $orphaned_2\nPaired: $paired\n";
```





### CNVnator_merger.py

```python
#!/usr/bin/env python3.5

'''
___________________________________________________

I take two CNVnator files called with different bins (e.g.: 500 bp and 1000 bp) and check concordance, significance and overlap and I print a unified file.
___________________________________________________

Andrea Del Cortona
2019/10/28
'''



#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import functools
import os
import re
import shutil
import string
import sys
import itertools
import operator
from datetime import datetime



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''I take two CNVnator files called with different bins (e.g.: 500 bp and 1000 bp) and check concordance, significance and overlap and I print a unified file.

The input files look like:
deletion	x1156_PM_chrIV:2001-17000	15000	0.00848363	1.06248e-11	0	1.22594e-11	0	1
deletion	x1156_PM_chrIV:44001-82000	38000	0.0501951	4.19401e-12	0	4.42701e-12	0	1

	usage:
	python3.5 CNVnator_merger.py --input_1 CNVnator_bin1 --input_2 CNVnator_bin2 --sample sample_name > CNV_merged.tab''')

parser.add_argument("--input_1",
	metavar ='INPUT_1',
	action = 'store',
	type = str,
	dest = 'INPUT_1',
	help = 'CNVnator output with bin1.',
	required = True)
	
parser.add_argument("--input_2",
	metavar ='INPUT_2',
	action = 'store',
	type = str,
	dest = 'INPUT_2',
	help = 'CNVnator output with bin2.',
	required = True)
	
parser.add_argument("--sample",
	metavar ='SAMPLE',
	action = 'store',
	type = str,
	dest = 'SAMPLE',
	help = 'The name of the sample.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

# read the input table and print the output
def main():

	# load databases
	with open(args.INPUT_1) as infile:
 		INPUT_1_DB = read_table(infile)

	with open(args.INPUT_2) as infile:
 		INPUT_2_DB = read_table(infile)

	# check positions, chromosome by chromosome. Print overlaps
	for KEY in INPUT_1_DB:
		if KEY in INPUT_2_DB:
			POS_1 = []
			POS_2 = []
			for k in range(len(INPUT_1_DB[KEY])):
				POS_1.append(INPUT_1_DB[KEY][k][0:2])									# here I save the positions of each record
			for k in range(len(INPUT_2_DB[KEY])):
				POS_2.append(INPUT_2_DB[KEY][k][0:2])									# here I save the positions of each record
			
			# iterate through positions and find overlapping regions
			for k in range(len(POS_1)):													# window 1
				for j in range(len(POS_2)):												# window 2
					if int(POS_1[k][1]) < int(POS_2[j][0]):								# window 1 < window 2
						continue
					elif int(POS_1[k][0]) > int(POS_2[j][1]):							# window 1 > window 2
						continue
					else:																# the windows overlap
						if int(POS_1[k][0]) <= int(POS_2[j][0]):
							if int(POS_1[k][1]) <= int(POS_2[j][1]):
								virtual_printer(args.SAMPLE, KEY, POS_2[j][0], POS_1[k][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
							else:
								virtual_printer(args.SAMPLE, KEY, POS_2[j][0], POS_2[j][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
						else:
							if int(POS_1[k][1]) <= int(POS_2[j][1]):
								virtual_printer(args.SAMPLE, KEY, POS_1[k][0], POS_1[k][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])
							else:
								virtual_printer(args.SAMPLE, KEY, POS_1[k][0], POS_2[j][1], INPUT_1_DB[KEY][k][2], INPUT_1_DB[KEY][k][3], INPUT_2_DB[KEY][j][3])



# import the input files
def read_table(infile):			

	IN_DB = {}

	for line in infile:
		line = line.rstrip('\n')
		KIND, LOCUS, LENGTH, NORM_RD, EVAL1, EVAL2, EVAL3, EVAL4, QUAL = line.split('\t')
		CHR, COHORD = LOCUS.split(':')
		START, STOP = COHORD.split('-')
		
		if float(EVAL1) < 0.05:
			if CHR in IN_DB:
				IN_DB[CHR].append([START, STOP, KIND, NORM_RD, EVAL1])
			else:
				IN_DB[CHR] = [[START, STOP, KIND, NORM_RD, EVAL1]]
	
	return(IN_DB)
	
	

# print the output
def virtual_printer(SAMPLE, CHR, START, STOP, KIND, NORM_RD1, NORM_RD2):
	MEAN_RD = float((float(NORM_RD1) + float(NORM_RD2))/2)
	print('\t'.join([SAMPLE, CHR, START, STOP, KIND, str(MEAN_RD)]))
	return()
	
	
	
#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
```







### Ping_MaltoseGenesParser.py

```python
#!/usr/bin/env python3.5

'''
___________________________________________________

I take the output of combined CNVnator and snpEff Annotated vcf files
contaning the 103 MAL negative strains of the study She et al.
https://doi.org/10.1016/j.cell.2017.12.015
and I print a strain-wise collection of info to understand the possible
common signature of loss of MAL genes.
___________________________________________________

Andrea Del Cortona
2019/11/28
'''


#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import functools
import os
import re
import shutil
import string
import sys
import itertools
import operator
import statistics
from scipy import stats
from statistics import mean
from datetime import datetime
from math import log10
from statsmodels.stats.multitest import multipletests



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''
#######################################################################

I take the output of combined CNVnator and snpEff Annotated vcf files
contaning the 103 MAL negative strains of the study She et al.
https://doi.org/10.1016/j.cell.2017.12.015
and I print a strain-wise collection of info to understand the possible
common signature of loss of MAL genes.

--> MAL_neg_strains.lst
	"3"	"SRR5629783"
	"7"	"SRR5629787"

--> S288C.MAltose_genes.lst
	ref_NC_001134_	800523	801929	YBR297W	MAL31
	ref_NC_001134_	802631	804475	YBR298C	MAL32

--> S288C.allMALneg_samples.tab
	S288C.SRR5629783	ref_NC_001133_	3001	24000	deletion	0.10224905000000001
	S288C.SRR5629783	ref_NC_001133_	27001	160000	duplication	1.43677

#######################################################################

	usage:
	python3.6 Ping_MaltoseGenesParser.py --strainlist MAL_neg_strains.lst \
		--Rlist RM11.Maltose_genes.lst \
		--Slist S288C.Maltose_genes.lst \
		--Ylist YJM975.Maltose_genes.lst \
		--Rvar 03_VARIANTS_RM11/RM11.combinedgVCFs.gen.MAL.snpEff.vcf \
		--Svar 03_VARIANTS_S288C/S288C.combinedgVCFs.gen.MAL.snpEff.vcf \
		--Yvar 03_VARIANTS_YJM975/YJM975.combinedgVCFs.gen.MAL.snpEff.vcf \
		--Rcnv 10_CNVs_RM11/RM11.allMALneg_samples.tab \
		--Scnv 10_CNVs_S288C/S288C.allMALneg_samples.tab \
		--Ycnv 10_CNVs_YJM975/YJM975.allMALneg_samples.tab
	
#######################################################################
	''')

parser.add_argument("--strainlist",
	metavar ='STRLIST',
	action = 'store',
	type = str,
	dest = 'STRLIST',
	help = 'SRR to numbered code list, aka MAL_neg_strains.lst.',
	required = True)
	
parser.add_argument("--Rlist",
	metavar ='RLIST',
	action = 'store',
	type = str,
	dest = 'RLIST',
	help = 'File with position of MAL genes in RM11, aka RM11.Maltose_genes.lst.',
	required = True)
	
parser.add_argument("--Slist",
	metavar ='SLIST',
	action = 'store',
	type = str,
	dest = 'SLIST',
	help = 'File with position of MAL genes in S288C, aka S288C.Maltose_genes.lst.',
	required = True)
			
parser.add_argument("--Ylist",
	metavar ='YLIST',
	action = 'store',
	type = str,
	dest = 'YLIST',
	help = 'File with position of MAL genes in YJM975, aka YJM975.Maltose_genes.lst.',
	required = True)

parser.add_argument("--Rvar",
	metavar ='RVAR',
	action = 'store',
	type = str,
	dest = 'RVAR',
	help = 'SnpEff.vcf file with annotated variants in MAL genes in RM11, aka RM11.combinedgVCFs.gen.MAL.snpEff.vcf.',
	required = True)
	
parser.add_argument("--Svar",
	metavar ='SVAR',
	action = 'store',
	type = str,
	dest = 'SVAR',
	help = 'SnpEff.vcf file with annotated variants in MAL genes in S288C, aka S288C.combinedgVCFs.gen.MAL.snpEff.vcf.',
	required = True)
	
parser.add_argument("--Yvar",
	metavar ='YVAR',
	action = 'store',
	type = str,
	dest = 'YVAR',
	help = 'SnpEff.vcf file with annotated variants in MAL genes in YJM975, aka YJM975.combinedgVCFs.gen.MAL.snpEff.vcf.',
	required = True)
	
parser.add_argument("--Rcnv",
	metavar ='RCNV',
	action = 'store',
	type = str,
	dest = 'RCNV',
	help = 'The CNVnator collapsed file, aka RM11.allMALneg_samples.tab.',
	required = True)
	
parser.add_argument("--Scnv",
	metavar ='SCNV',
	action = 'store',
	type = str,
	dest = 'SCNV',
	help = 'The CNVnator collapsed file, aka S288C.allMALneg_samples.tab.',
	required = True)
	
parser.add_argument("--Ycnv",
	metavar ='YCNV',
	action = 'store',
	type = str,
	dest = 'YCNV',
	help = 'The CNVnator collapsed file, aka YJM975.allMALneg_samples.tab.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

#------------------------------------------------------------------#
# read the input tables and print the output
def main():
	
	# import strain list
	STRAIN_DB = []
	with open(args.STRLIST) as listfile:
		for line in listfile:
			line = line.strip('\n')
			STRNUM, STRAIN = line.split('\t')
			STRAIN_DB.append([STRNUM.replace('"', ''), STRAIN.replace('"', '')])

	# import MAL genes positions
	RPOS_DB = import_positions(args.RLIST)
	SPOS_DB = import_positions(args.SLIST)
	YPOS_DB = import_positions(args.YLIST)
	
	# import snpEff.vcf files
	RVAR_DB = import_variants(args.RVAR)
	SVAR_DB = import_variants(args.SVAR)
	YVAR_DB = import_variants(args.YVAR)
	
	# import merged CNVs files
	RCNV_DB = import_cnvs(args.RCNV)
	SCNV_DB = import_cnvs(args.SCNV)
	YCNV_DB = import_cnvs(args.YCNV)
	
	# process samples
	RPROC_DB = processdb(STRAIN_DB, RPOS_DB, RVAR_DB, RCNV_DB)
	SPROC_DB = processdb(STRAIN_DB, SPOS_DB, SVAR_DB, SCNV_DB)
	YPROC_DB = processdb(STRAIN_DB, YPOS_DB, YVAR_DB, YCNV_DB)

	# print output
	GENELIST = ["IMA1", "IMA2", "IMA3", "IMA4", "IMA5", "MAL11", "MAL12", "MAL13", "MAL31", "MAL32", "MAL33", "Maltose_responsive_factor", "MPH2", "MPH3"]
	HEADER = [[], [], [], [], [], []]
	OUTPUT = []

	for SAMPLE in STRAIN_DB:
		OUTLINE = []	
					
		# first print variants
		for GENE in GENELIST:
			if GENE not in HEADER[0]:
				HEADER[0].append(GENE)
			for ELEMENT in RPROC_DB[0][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))
		for GENE in GENELIST:
			if GENE not in HEADER[1]:
				HEADER[1].append(GENE)
			for ELEMENT in SPROC_DB[0][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))
		for GENE in GENELIST:				
			if GENE not in HEADER[2]:
				HEADER[2].append(GENE)
			for ELEMENT in YPROC_DB[0][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))					

		# now print CNVs
		for GENE in GENELIST:
			if GENE not in HEADER[3]:
				HEADER[3].append(GENE)
			for ELEMENT in RPROC_DB[1][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))
		for GENE in GENELIST:
			if GENE not in HEADER[4]:
				HEADER[4].append(GENE)
			for ELEMENT in SPROC_DB[1][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))
		for GENE in GENELIST:
			if GENE not in HEADER[5]:
				HEADER[5].append(GENE)
			for ELEMENT in YPROC_DB[1][SAMPLE[0], SAMPLE[1]]:
				ALMOST = []
				if GENE == ELEMENT[0]:
					for ENTRY in ELEMENT[1]:
						ALMOST.append(','.join(str(X) for X in ENTRY))
					OUTLINE.append(" ".join(ALMOST))

		OUTPUT.append("\t".join([SAMPLE[0], SAMPLE[1], "\t".join(OUTLINE)]))
			
	A = "\t".join(HEADER[0])
	B = "\t".join(HEADER[1])
	C = "\t".join(HEADER[2])
	D = "\t".join(HEADER[3])
	E = "\t".join(HEADER[4])
	F = "\t".join(HEADER[5])
	
	print('\t'.join(["# Sample_number", "Sample_SRR", A, B, C, D, E, F]))
	for LINE in OUTPUT:
		print(LINE)
#------------------------------------------------------------------#



#------------------------------------------------------------------#
# I import the MAL position files
def import_positions(LIST):
	
	OUTLIST = []
	with open(LIST) as listfile:
		for line in listfile:
			line = line.strip('\n')
			CHR, START, STOP, SGDCODE, GENENAME = line.split('\t')
			OUTLIST.append([CHR, START, STOP, SGDCODE, GENENAME])
	
	return(OUTLIST)
#------------------------------------------------------------------#



#------------------------------------------------------------------#
# I import the annotated variants
def import_variants(LIST):
	
	OUTLIST = []
	with open(LIST) as listfile:
		for line in listfile:
			line = line.strip('\n')
			if re.match('#', line):
				continue
			else:
				CHR, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *SAMPLES = line.split('\t')
				OUTLIST.append([CHR, POS, REF, ALT, INFO, FORMAT, SAMPLES])
	
	return(OUTLIST)
#------------------------------------------------------------------#



#------------------------------------------------------------------#
# I import the merged CNVnator files
def import_cnvs(LIST):
	
	OUTLIST = []
	with open(LIST) as listfile:
		for line in listfile:
			line = line.strip('\n')
			RF_SAMPLE, CHR, START, STOP, KIND, QUANT = line.split('\t')
			REF, SAMPLE = RF_SAMPLE.split('.')
			OUTLIST.append([SAMPLE, CHR, START, STOP, KIND, QUANT])
	
	return(OUTLIST)
#------------------------------------------------------------------#



#------------------------------------------------------------------#
# I process the samples
def processdb(STRAIN_DB, POS_DB, VAR_DB, CNV_DB):
	
	# initiliaze variables
	GENELIST = ["IMA1", "IMA2", "IMA3", "IMA4", "IMA5", "MAL11", "MAL12", "MAL13", "MAL31", "MAL32", "MAL33", "Maltose_responsive_factor", "MPH2", "MPH3"]
	OUTVAR = {}
	OUTCNV = {}
	k = 0
	
	# iterate through samples
	for SAMPLE in STRAIN_DB:
		STRNUM = SAMPLE[0]
		STRAIN = SAMPLE[1]
		VAR_OUT = {}
		CNV_OUT = {}

		# iterate through positions
		for POS_ENTRY in POS_DB:
			CHR = POS_ENTRY[0]
			PSTART = POS_ENTRY[1]
			PSTOP = POS_ENTRY[2]
			SGDCODE = POS_ENTRY[3]
			GENENAME = POS_ENTRY[4]
			
			# integrate variants
			for VAR_ENTRY in VAR_DB:
				PCHR = VAR_ENTRY[0]
				POS = VAR_ENTRY[1]
				REF = VAR_ENTRY[2]
				ALT = VAR_ENTRY[3]
				INFO = VAR_ENTRY[4]
				FORMAT = VAR_ENTRY[5]
				SAMPLES = VAR_ENTRY[6]
			
				if CHR == PCHR and int(POS) > int(PSTART) and int(POS) < int(PSTOP):
					ANNOT, KIND, EFF, *MORE = INFO.split('|')
					GT, *ANN = SAMPLES[k].split(':')
					GT_TEST = 0
					
					# genotype check
					if GT == "0/0" or GT == "0|0":
						GT_TEST = 0
					else:
						GT_TEST = 1
	
					# annotate variant for the sample
					if GENENAME in VAR_OUT:
						VAR_OUT[GENENAME].append([":".join([CHR, "-".join([PSTART, PSTOP])]), POS, ALT, ":".join([KIND, EFF]), GT_TEST])
					else:
						VAR_OUT[GENENAME] = [[":".join([CHR, "-".join([PSTART, PSTOP])]), POS, ALT, ":".join([KIND, EFF]), GT_TEST]]
			
			# integrate cnvs	
			for CNV_ENTRY in CNV_DB:
				if STRAIN in CNV_ENTRY:
					CSAMPLE = CNV_ENTRY[0]
					CCHR = CNV_ENTRY[1]
					CSTART = CNV_ENTRY[2]
					CSTOP = CNV_ENTRY[3]
					CKIND = CNV_ENTRY[4]
					CQUANT = CNV_ENTRY[5]
					
					# same chromosome
					if CCHR == CHR:
					
						# do windows overlap?
						if CSTOP < PSTART:
							continue
						elif CSTART > PSTOP:
							continue
						else:

							# annotate cnv for the sample
							if GENENAME in CNV_OUT:
								CNV_OUT[GENENAME].append([":".join([CHR, "-".join([PSTART, PSTOP])]), ":".join([CKIND, CQUANT])])
							else:
								CNV_OUT[GENENAME] = [[":".join([CHR, "-".join([PSTART, PSTOP])]), ":".join([CKIND, CQUANT])]]

		# prepare output
		OUTVAR[SAMPLE[0], SAMPLE[1]] = []
		OUTCNV[SAMPLE[0], SAMPLE[1]] = []
		for GENE in GENELIST:
			if GENE in VAR_OUT:
				OUTVAR[SAMPLE[0], SAMPLE[1]].append([GENE, VAR_OUT[GENE]])
			else:
				OUTVAR[SAMPLE[0], SAMPLE[1]].append([GENE, ''])
			
			if GENE in CNV_OUT:
				OUTCNV[SAMPLE[0], SAMPLE[1]].append([GENE, CNV_OUT[GENE]])
			else:
				OUTCNV[SAMPLE[0], SAMPLE[1]].append([GENE, ''])
				
		k += 1	
		
	return(OUTVAR, OUTCNV)
#------------------------------------------------------------------#

	

#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
```





### Ping_IMA1_SNPs_simplified.py

```python
#!/usr/bin/env python3.5

'''
___________________________________________________

I take a genotyped vcf file with the positions around IMA1 locus
and I transform it into a manageable table for plotting
___________________________________________________

Andrea Del Cortona
2020/03/02
'''



#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import functools
import os
import re
import shutil
import string
import sys
import itertools
import operator
import statistics
from scipy import stats
from statistics import mean
from datetime import datetime
from math import log10
from statsmodels.stats.multitest import multipletests



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''
#######################################################################

I take a genotyped vcf file with the positions around IMA1 locus
and I transform it into a manageable table for plotting

#######################################################################

	usage:
	python3.6 Ping_IMA1_SNPs_simplified.py --input combinedGVCFs.vcf > combinedGVCFs.tab
	
#######################################################################
	''')

parser.add_argument("--input",
	metavar ='INFILE',
	action = 'store',
	type = str,
	dest = 'INFILE',
	help = 'genotyped VCF file with IMA1 positions.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

#------------------------------------------------------------------#
# read the input tables and print the output
def main():
	
	# open input file
	POS_LIST = ["SAMPLE"]
	SAMPLE_LIST = []
	SAMPLE_DB = {}
	
	with open(args.INFILE) as infile:
		for line in infile:
			line = line.strip('\n')
			if re.match("##", line):
				# skip header
				continue
				
			elif re.match("#", line):
				# save sample names
				CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *MORE = line.split('\t')
				for SAMPLE in MORE:
					SAMPLE_LIST.append(SAMPLE)
					
			else:
				# process position
				# we have two FORMAT styles:
				# GT:AD:DP:GQ:PGT:PID:PL and GT:AD:DP:GQ:PL
				# I need to take this into account
				#
				# I have to take into account as well positions with no coverage
				# GT = ./.				
				CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *MORE = line.split('\t')	
				POS_LIST.append(':'.join([CHROM, POS]))
				
				for k in range(len(SAMPLE_LIST)):
				
					GT, AD, DP, GQ, *EVENMORE = MORE[k].split(':')
					if GT == "0/0":
						ALT1 = REF
					elif GT == "./.":
						ALT1 = "NO_COV"
					else:
						ALT1 = ALT
					
					if SAMPLE_LIST[k] in SAMPLE_DB:
						SAMPLE_DB[SAMPLE_LIST[k]].append(ALT1)
					else:
						SAMPLE_DB[SAMPLE_LIST[k]] = [ALT1]
				
	print('\t'.join(POS_LIST))
	for key in SAMPLE_DB:
		print('\t'.join([key, '\t'.join(SAMPLE_DB[key])]))
#------------------------------------------------------------------#

	

#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
```





### Ping_IMA1_CH408044.30kb_cov.py

```python
#!/usr/bin/env python3.5

'''
___________________________________________________

I take a the coverge info for the first 30kb of chromosome CH408044
and I calculate the average coverage in the window 10-20 kb.
I identify samples with coverage close to 0 in this region.
___________________________________________________

Andrea Del Cortona
2020/04/01
'''



#==================================================================#
#   LOAD LIBRARIES                                                 #
#==================================================================#

import argparse
import functools
import os
import re
import shutil
import string
import sys
import itertools
import operator
import statistics
from scipy import stats
from statistics import mean
from datetime import datetime
from math import log10
from statsmodels.stats.multitest import multipletests



#==================================================================#
#   INPUT PARSER                                                   #
#==================================================================#

parser = argparse.ArgumentParser(description='''
#######################################################################

I take a the coverge info for the first 30kb of chromosome CH408044
and I calculate the average coverage in the window 10-20 kb.
I identify samples with coverage close to 0 in this region.

#######################################################################

	usage:
	python3.6 Ping_IMA1_CH408044.30kb_cov.py --input /media/DISK2-3TB/Ping/strain_file_CORRECT20200424.txt
	
#######################################################################
	''')

parser.add_argument("--input",
	metavar ='INFILE',
	action = 'store',
	type = str,
	dest = 'INFILE',
	help = '/media/DISK2-3TB/Ping/strain_file_CORRECT20200424.txt strain list file. On the first columns it has the strain number, on the second column the corresponding SRR number.',
	required = True)

args = parser.parse_args()



#==================================================================#
#   FUNCTIONS                                                      #
#==================================================================#

#------------------------------------------------------------------#
# read the input tables and print the output
def main():
	
	# read strain list
	STRAIN_DB = []
	with open(args.INFILE) as infile:
		for line in infile:
			line = line.strip('\n')
			NUM, SRR = line.split('\t')
			STRAIN_DB.append([NUM.replace('"', ''), SRR.replace('"', '')])
	
	# open coverage files
	OUT_FILE = []
	for k in range(len(STRAIN_DB)):
		NUM = STRAIN_DB[k][0]
		SRR = STRAIN_DB[k][1]
		NAME_FILE =  SRR+".CH408044.30kb.1kb_cov.bed"
		
		with open(NAME_FILE) as covfile:
			# read the coverage
			FILE_POS = []
			for line in covfile:
				line = line.strip('\n')
				CHR, START, STOP, COV = line.split('\t')
				FILE_POS.append([CHR, START, STOP, COV])	
			
			# select coverage in win 10-20 kb and compute average
			COV_WIN = []
			for i in range(10,20):
				COV_WIN.append(float(FILE_POS[i][3]))
			COV_AVG = mean(COV_WIN)
			if COV_AVG < 1:
				FLAG = "SHORT"
			else:
				FLAG = "LONG"
		
		OUT_FILE.append([SRR, NUM, str(round(COV_AVG, 4)), FLAG])
		
	# print output
	for k in OUT_FILE:
		print('\t'.join(k))	
#------------------------------------------------------------------#

	

#==================================================================#
#   RUN                                                            #
#==================================================================#

# run script and give the running time
if __name__ == '__main__':
	t0 = datetime.now()
	main()
	dt = datetime.now() - t0
sys.stderr.write( "# Time elapsed: %s\n" % dt )
```





### MALgene_combinedmatrix.R

```R
## load libraries
library(ggplot2)
library(gridExtra)
library(plyr)

# upload files
setwd("/media/DISK2-3TB/Ping/30_S288C_IMA1_allele_correct20200424/")
MALgenes = read.delim("strains_MAL_gen2phe.txt", header = TRUE)

#------------------------------------------------------------------#
# IMA1 allele length
MALgenes$IMA1_ALLELE = factor(MALgenes$IMA1_ALLELE, levels=c("LONG", "SHORT", "FAILED"))
mu = ddply(drop_na(MALgenes), "IMA1_ALLELE", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = IMA1_ALLELE, color = IMA1_ALLELE)) +
  geom_histogram(bins = 20, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  ylim(0, 120) +
  coord_cartesian(expand = FALSE) +
  labs(title = "IMA1 allele and Maltose growth",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 28, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = c(0.9, 0.85),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# IMA1
# create output grid layout
grid_layout <- rbind(c(1, 2),
                     c(3, 4))

# SNV_4835
MALgenes$SNV4835 = factor(MALgenes$SNV4835, levels=c("G", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV4835", summarise, grp.mean=mean(Maltose_growth))
SNV_4835 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV4835, color = SNV4835)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV4835",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_4836
MALgenes$SNV4836 = factor(MALgenes$SNV4836, levels=c("T", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV4836", summarise, grp.mean=mean(Maltose_growth))
SNV_4836 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV4836, color = SNV4836)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV4836",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_4837
MALgenes$SNV4837 = factor(MALgenes$SNV4837, levels=c("A", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV4837", summarise, grp.mean=mean(Maltose_growth))
SNV_4837 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV4837, color = SNV4837)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV4837",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_4838
MALgenes$SNV4838 = factor(MALgenes$SNV4838, levels=c("C", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV4838", summarise, grp.mean=mean(Maltose_growth))
SNV_4838 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV4838, color = SNV4838)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV4838",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# print the output graphs
grid.arrange(SNV_4835, SNV_4836, SNV_4837, SNV_4838,
             layout_matrix = grid_layout,
             top =textGrob("IMA1", gp=gpar(fontsize=24)))

#------------------------------------------------------------------#
# IMA5
# SNV5983
MALgenes$SNV5983 = factor(MALgenes$SNV5983, levels=c("C", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV5983", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV5983, color = SNV5983)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "IMA5: SNV5983",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# MAL13
# VII:1070598
MALgenes$VII.1070598 = factor(MALgenes$VII.1070598, levels=c("A", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "VII.1070598", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = VII.1070598, color = VII.1070598)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "MAL13: VII:1070598",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# MAL11
# VII:1074352
MALgenes$VII.1074352 = factor(MALgenes$VII.1074352, levels=c("T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "VII.1074352", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = VII.1074352, color = VII.1074352)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "MAL11: VII:1074352",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# MAL31
#create output grid layout
grid_layout <- rbind(c(1, 2),
                     c(3, 4))

# SNV_884
MALgenes$SNV884 = factor(MALgenes$SNV884, levels=c("T", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV884", summarise, grp.mean=mean(Maltose_growth))
SNV_884 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV884, color = SNV884)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV884",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_885
MALgenes$SNV885 = factor(MALgenes$SNV885, levels=c("T", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV885", summarise, grp.mean=mean(Maltose_growth))
SNV_885 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV885, color = SNV885)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV885",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_886
MALgenes$SNV886 = factor(MALgenes$SNV886, levels=c("G", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV886", summarise, grp.mean=mean(Maltose_growth))
SNV_886 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV886, color = SNV886)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV886",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# print the output graphs
blank <- grid.rect(gp=gpar(col="white"))
grid.arrange(SNV_884, SNV_885, SNV_886, blank,
             layout_matrix = grid_layout,
             top =textGrob("MAL31", gp=gpar(fontsize=24)))

#------------------------------------------------------------------#
# MAL32
# II:802702
MALgenes$II.802702 = factor(MALgenes$II.802702, levels=c("T", "C", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "II.802702", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = II.802702, color = II.802702)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "MAL32: II:802702",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# MAL33
#create output grid layout
grid_layout <- rbind(c(1, 2),
                     c(3, 4))

# SNV_872
MALgenes$SNV872 = factor(MALgenes$SNV872, levels=c("A", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV872", summarise, grp.mean=mean(Maltose_growth))
SNV_872 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV872, color = SNV872)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV872",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_873
MALgenes$SNV873= factor(MALgenes$SNV873, levels=c("T", "A", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV873", summarise, grp.mean=mean(Maltose_growth))
SNV_873 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV873, color = SNV873)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV873",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_875
MALgenes$SNV875 = factor(MALgenes$SNV875, levels=c("G", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV875", summarise, grp.mean=mean(Maltose_growth))
SNV_875 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV875, color = SNV875)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV875",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# SNV_876
MALgenes$SNV876 = factor(MALgenes$SNV876, levels=c("A", "G", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV876", summarise, grp.mean=mean(Maltose_growth))
SNV_876 = ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV876, color = SNV876)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "SNV876",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# print the output graphs
grid.arrange(SNV_872, SNV_873, SNV_875, SNV_876,
             layout_matrix = grid_layout,
             top =textGrob("MAL33", gp=gpar(fontsize=24)))

#------------------------------------------------------------------#
# MPH2
# IV:7726
MALgenes$IV.7726 = factor(MALgenes$IV.7726, levels=c("C", "A", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "IV.7726", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = IV.7726, color = IV.7726)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "MPH2: IV:7726",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# MPH3
# X:738075
MALgenes$X.738075= factor(MALgenes$X.738075, levels=c("T", "G", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "X.738075", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = X.738075, color = X.738075)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "MPH3: X:738075",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

#------------------------------------------------------------------#
# YPR196W
# SNV11999
MALgenes$SNV11999 = factor(MALgenes$SNV11999, levels=c("A", "T", "FAILED", "NO_COV"))
mu = ddply(drop_na(MALgenes), "SNV11999", summarise, grp.mean=mean(Maltose_growth))
ggplot(MALgenes, aes(x=Maltose_growth, fill = SNV11999, color = SNV11999)) +
  geom_histogram(bins = 25, position="identity", alpha=0.2) +
  geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed", colour=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "pink")) +
  ylim(0, 100) +
  coord_cartesian(expand = FALSE) +
  labs(title = "YPR196W: SNV11999",
       x = "Maltose growth",
       y = "Number of segregants") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title = element_text(size = 18),
        axis.text= element_text(size = 14),
        legend.background = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
```





# _____________________________________________________________________________________________



## - Andrea Del Cortona -

## - 22 September 2020 -



# _____________________________________________________________________________________________

## 

