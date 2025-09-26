# recombination-landscape
Preparing VCF and running ReLERNN

### 1. For each indiviudal, run haplotypecaller using GATK to produce a GVCF file

```
#for gpu compatible devices I use Parabricks

pbrun fq2bam \
  --ref Leiomod.fa \
  --in-fq Leiomadagascarensis-ran38059_S10_L004_R1_001.fastq.gz Leiomadagascarensis-ran38059_S10_L004_R2_001.fastq.gz \
  --out-bam LM_ran38059.bam

pbrun haplotypecaller \
  --ref Leiomod.fa \
  --in-bam LM_ran38059.bam \
  --out-variants LM_ran38059.g.vcf --gvcf
```

#### 2. Combine all with CombinGVCFs. Make sure the name of the sample is in the GVCFs

```
while IFS= read -r name; do
  # Your command here, using $line for the current line
  echo "Processing line: $name"
  bgzip ${name}.g.vcf
  tabix ${name}.g.vcf.gz
  echo ${name} > ${name}
  bcftools reheader -s $name -o ${name}.rename.g.vcf.gz ${name}.g.vcf.gz
  tabix -f ${name}.rename.g.vcf.gz
done < leionames.txt

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /mendel-nas1/ddebaun/gatk/build/libs/gatk-package-4.2.5.0-17-gb14f810-SNAPSHOT-local.jar CombineGVCFs -R $reference -V LLM_RAN39039.rename.g.vcf.gz -V LLM_RAX12322.rename.g.vcf.gz -V LLM_RAX12546.rename.g.vcf.gz -V LLM_RAX5299.rename.g.vcf.gz -V LLM_ran38059.rename.g.vcf.gz -V LLM_rax12515.rename.g.vcf.gz -O Leiomod_combined.g.vcf.gz

```
#### 3. Run GenotypeGVCF 

```
java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /mendel-nas1/ddebaun/gatk/build/libs/gatk-package-4.2.5.0-17-gb14f810-SNAPSHOT-local.jar GenotypeGVCFs -R $reference -V Leiomod_combined.g.vcf.gz -O Leiomod_allsites.vcf.gz --include-non-variant-sites

```
#### 4. Filter VCF file

First, get all biallelic SNPs
```
vcf_in="${name}.vcf"
filtered_vcf="${name}_snps.vcf"
bcftools view \
  -f .,PASS \
  -v snps \
  --min-alleles 2 --max-alleles 2 \
  $vcf_in \
  -Ov -o $filtered_vcf
```

For samples with low depth and low genotype quality, convert these to missing.
```
masked_vcf="${name}_masked.vcf.gz"
bcftools +setGT $filtered_vcf -- -t q -n . -i 'FMT/DP < 5 || FMT/GQ < 30' \
| bcftools view -Oz -o $masked_vcf
```

ReLERNN learns from the pattern of polymorphism across haplotypes. Very few samples (or sites present in only 1–2 samples) give weak signal and increase noise. Here I require no more than 33% of samples are missing. For 6 taxa, requires 4 are present at a site.
```
filtered_vcf_2="${name}_snps_relernn.vcf.gz"
bcftools view \
       -i "F_MISSING <= 0.33" \
       $masked_vcf \
        -Oz -o $filtered_vcf_2
```

Additional quality control filtering of SNPs 
```
tabix ${filtered_vcf_2}

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /mendel-nas1/ddebaun/gatk/build/libs/gatk-package-4.2.5.0-17-gb14f810-SNAPSHOT-local.jar VariantFiltration \
-R $reference \
-V ${filtered_vcf_2}.gz \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "ExcessHet > 54.69" --filter-name "ExcessHet54.69" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O ${name}_filtered.vcf.gz

bcftools view  -f 'PASS' ${name}_filtered.vcf.gz  -Ov -o ${name}_filtered_PASS.vcf
```

#### 5. Create mask file
Some sites of the reference genome may never be mapped to across your samples. Here, I create a mask of these regions based on depth and genotype quality statistics

```
vcftools --vcf Leiomod_allsites.vcf.gz --site-mean-depth --out Leiomod_depth
vcftools --vcf Leiomod_allsites.vcf.gz --site-quality --out Leiomod

awk '$4 < 5 {print $1"\t"$2-1"\t"$2}' Leiomod_depth.ldepth.mean > Leiomod_lowdepth.bed
awk '$3 < 30 {print $1"\t"$2-1"\t"$2}' Leiomod.lqual > Leiomod_lowqual.bed

cat Leiomod_lowdepth.bed Leiomod_lowqual.bed | sort -k 1,1 -k2,2n | bedtools merge > Leiomod.mask
```


#### 6. Create genome file
```
samtools faidx ../Leiomod.fa
cut -f1,2 ../Leiomod.fa.fai > Leiomod.chrom.sizes
awk '{print $1"\t0\t"$2}' Leiomod.chrom.sizes > Leiomod.chrom.bed

```

### 7. Run ReLERNN

Simulating based on the VCF
```
name="Leiomod"
DIR=${name}_results
mkdir $DIR
VCF=${name}_ReLERNN.vcf
GENOME=Leiomod.chrom.bed
MU="2.4e-9" #snake mu from green et al.
ReLERNN_SIMULATE --vcf ${VCF} --genome ${GENOME} --unphased --projectDir ${DIR} --assumedMu ${MU} --assumedGenTime 2 --nCPU $SLURM_NTASKS_PER_NODE --mask ${name}.mask
```

Training, Prediction, and Bootstrapping
```
ReLERNN_TRAIN --projectDir ${DIR}
ReLERNN_PREDICT --vcf ${VCF} --projectDir ${DIR} --minSites 150 --unphased 
ReLERNN_BSCORRECT --projectDir ${DIR} --nCPU $SLURM_NTASKS_PER_NODE --nSlice 100 --nReps 1000
```








