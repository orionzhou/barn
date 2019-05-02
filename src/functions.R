require(devtools)
load_all('~/git/rmaize')
require(rentrez)
require(xml2)
dird = '~/projects/barn/data'
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
    if('directory' %in% colnames(ti)) ti = ti %>% fill(directory)
    ti %>% fill(Tissue, Genotype, Treatment) %>%
        select(SampleID, Tissue, Genotype, Treatment, Replicate, everything()) %>%
        arrange(SampleID, Tissue, Genotype, Treatment) %>%
        group_by(Tissue,Genotype,Treatment) %>%
        mutate(Replicate = 1:n()) %>%
        ungroup()
    #}}}
}

locate_fastq <- function(diri, file_prefix, fmt='umgc', interleaved=F) {
    #{{{
    if(fmt == 'umgc') {
        r0 = sprintf("%s_R1_001.fastq", file_prefix)
        r1 = sprintf("%s_R1_001.fastq", file_prefix)
        r2 = sprintf("%s_R2_001.fastq", file_prefix)
    } else if (fmt == 'jgi') {
        r0 = sprintf("%s.fastq", file_prefix)
        r1 = sprintf("%s_1.fastq", file_prefix)
        r2 = sprintf("%s_2.fastq", file_prefix)
    } else if (fmt == 'simple') {
        r0 = sprintf("%s.fq", file_prefix)
        r1 = sprintf("%s_1.fq", file_prefix)
        r2 = sprintf("%s_2.fq", file_prefix)
    } else {
        stop(sprintf("unknown fmt: %s\n", opt))
    }
    r0 = file.path(diri, r0)
    r1 = file.path(diri, r1)
    r2 = file.path(diri, r2)
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




