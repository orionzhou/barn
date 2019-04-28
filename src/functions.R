require(devtools)
load_all('~/git/rmaize')
dird = '~/projects/local/data'
f_cfg = file.path(dird, '01.cfg.xlsx')
t_cfg = read_xlsx(f_cfg, sheet=1, col_names=T) %>%
    mutate(interleaved = as.logical(interleaved)) %>%
    mutate(done = as.logical(done)) %>%
    mutate(run = as.logical(run)) %>%
    replace_na(list(interleaved=F, done=F, run=F))

generate_sample_id <- function(n, pre='s') {
    #{{{
    ndig = floor(log10(n)) + 1
    fmt = sprintf("s%%0%dd", ndig)
    sprintf(fmt, 1:n)
    #}}}
}

complete_sample_list <- function(ti) {
    #{{{
    if(!'Treatment' %in% colnames(ti)) ti = ti %>% mutate(Treatment='')
    if(!'Replicate' %in% colnames(ti)) ti = ti %>% mutate(Replicate='')
    if( is.na(ti$SampleID[1]) ) ti$SampleID = generate_sample_id(nrow(ti))
    ti %>% fill(Tissue, Genotype, Treatment, directory) %>%
        select(SampleID, Tissue, Genotype, Treatment, Replicate, everything())
    #}}}
}

locate_fastq <- function(diri, file_prefix, opt='umgc', interleaved=F) {
    #{{{
    if(opt == 'umgc') {
        r0 = sprintf("%s/%s_R1_001.fastq", diri, file_prefix)
        r1 = sprintf("%s/%s_R1_001.fastq", diri, file_prefix)
        r2 = sprintf("%s/%s_R2_001.fastq", diri, file_prefix)
    } else if (opt == 'jgi') {
        r0 = sprintf("%s/%s.fastq", diri, file_prefix)
        r1 = sprintf("%s/%s_1.fastq", diri, file_prefix)
        r2 = sprintf("%s/%s_2.fastq", diri, file_prefix)
    } else if (opt == '3rnaseq') {
        r0 = sprintf("%s/%s.fq", diri, file_prefix)
        r1 = sprintf("%s/%s_1.fq", diri, file_prefix)
        r2 = sprintf("%s/%s_2.fq", diri, file_prefix)
    } else {
        stop(sprintf("unknown opt: %s\n", opt))
    }
    if (!file.exists(r0) & !file.exists(r1) & !file.exists(r2)) {
        r0 = sprintf("%s.gz", r0)
        r1 = sprintf("%s.gz", r1)
        r2 = sprintf("%s.gz", r2)
    }
    if (file.exists(r1) & file.exists(r2)) {
        r0 = ''; paired = T
    } else if (file.exists(r0)) {
        r1 = ''; r2 = ''; paired = F
    } else {
        stop(sprintf("fastq not found for %s/%s\n", diri, file_prefix))
    }
    if (interleaved) {
        stopifnot(!paired)
        paired = T
    }
    tibble(paired=paired,interleaved=interleaved,r0=r0,r1=r1,r2=r2)
    #}}}
}

read_msi_fastqc <- function(rdir) {
    #{{{
    fi = file.path(rdir, "Analysis/illumina-basicQC/Resources/metrics.txt")
    if(file.exists(fi)) {
        read_tsv(fi) %>%
        select(sampleName = `general-samplename`,
               spots = `fastqc-totalsequences`,
               avgLength = `fastqc-sequencelength`) %>%
        filter(sampleName != 'Mean')
    } else{
        NA
    }
    #}}}
}



