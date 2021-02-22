A632 variants were called using the [nf-core/sarek pipeline](https://github.com/nf-core/sarek) (v2.7). In short, raw sequencing reads went through quality control with [FastQC], trimmed using [Trim-Galore] and mapped to the maize [B73 AGP_v4 genome] using [BWA-mem]. Duplicates were marked using [GATK MarkDuplicatesSpark].  Base quality scores were recalibrated using [GATK BaseRecalibrator] and [ApplyBQSR].  Raw variants were called using [GATK Haplotypcaller] and [GenotypeGVCF].  Final (high confidence) variants were obtained by applying a "homozygous" and "genotype quality" filter (i.e., `bcftools view -i 'GT="AA" && GQ>=20'`).

- [Pipeline report](https://s3.msi.umn.edu/zhoup-nfo/zm.dn21a/Reports/MultiQC/multiqc_report.html)
- [Raw VCF](https://s3.msi.umn.edu/zhoup-nfo/zm.dn21a/VariantCalling/A632/HaplotypeCaller/HaplotypeCaller_A632.vcf.gz) (including 12,389,787 SNPs + 1,429,519 InDels) and [index](https://s3.msi.umn.edu/zhoup-nfo/zm.dn21a/VariantCalling/A632/HaplotypeCaller/HaplotypeCaller_A632.vcf.gz.tbi)
- [Filtered (high quality) VCF](https://s3.msi.umn.edu/zhoup-share/zm.dn21a/A632.hq.vcf.gz) (5,587,291 SNPs + 864,901 InDels), [index](https://s3.msi.umn.edu/zhoup-share/zm.dn21a/A632.hq.vcf.gz.tbi) and [stats](https://s3.msi.umn.edu/zhoup-share/zm.dn21a/A632.hq.stats)











