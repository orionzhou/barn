source("functions.R")
t_cfg %>% count(libtype,source,format)
t_cfg %>% select(-stranded,-interleaved,-readtype) %>% print(n=50)
rcfg = t_cfg %>% filter(source=='sra')
lcfg = t_cfg %>% filter(source=='local')

#{{{ SRA - work on single one
yid = 'rn12a'
acc = rcfg %>% filter(yid == !!yid) %>% pull(accession)
ti = get_sra_meta(acc, yid)
fo = sprintf("%s/08_sra_list_raw/%s.csv", dird, yid)
write_csv(ti, fo)

#{{{ work on random project
acc = 'PRJNA305809'
ti = get_sra_meta(acc, '')
fo = sprintf("%s/tmp.csv", dird)
write_csv(ti, fo)
#}}}

fx = file.path(dird, '09.sra.xlsx')
tx = read_xlsx(fx, sheet=yid)
tx = complete_sample_list(tx)
fo = sprintf("%s/09_sra_list/%s.tsv", dird, yid)
write_tsv(tx, fo, na = '')
#}}}

#{{{ SRA- work on many
yids = rcfg %>% filter(!yid %in% c("rn18e",'rn99a')) %>% pull(yid)
yids = yids[yids %in% c('dn12a','dn17b')]
to = rcfg %>% select(yid, accession) %>%
    filter(yid %in% yids, !is.na(accession)) %>%
    mutate(fo = sprintf("%s/08_sra_list_raw/%s.csv", dird, yid)) %>%
    mutate(data = map2(accession, yid, get_sra_meta))
to %>% mutate(x = map2(data, fo, write_csv))

ti = to %>% filter(yid == !!yid) %>% select(data) %>% unnest()
ti %>% print(width=Inf)
ti %>% count(lib_layout,lib_selection,lib_source,lib_strategy)

# manually add to '09.sra.xlsx'
fx = file.path(dird, '09.sra.xlsx')
sheets = excel_sheets(fx)
tx = tibble(fx = fx, sheet = sheets) %>%
    mutate(fo = sprintf("%s/09_sra_list/%s.tsv", dird, sheet)) %>%
    mutate(ti = map2(fx, sheet, read_xlsx)) %>%
    mutate(to = map(ti, complete_sample_list))
tx %>% mutate(x = map2(to, fo, write_tsv, na = ''))
#}}}

#{{{ local data
yid = 'rn99f'
fmt = lcfg %>% filter(yid==!!yid) %>% pull(format)
lid = lcfg %>% filter(yid==!!yid) %>% pull(lid)
interleaved = lcfg %>% filter(yid==!!yid) %>% pull(interleaved)
fi = file.path(dird, '05.local.xlsx')
ti = read_xlsx(fi, sheet=lid)
to = complete_sample_list(ti) %>%
    mutate(data = map2(directory,file_prefix, locate_fastq,
                       fmt=!!fmt, interleaved=!!interleaved)) %>%
    unnest() %>% select(-directory,-file_prefix)
to %>% count(MergeID, Genotype, Tissue, Treatment) %>% print(n=50)
#
fo = sprintf("%s/06_local_list/%s.tsv", dird, yid)
write_tsv(to, fo)
#}}}

fix_read_list <- function(ti, yid) {
#{{{
if(yid == 'dn12a') {
    #{{{ hapmap2
    th = ti %>%
        transmute(SampleID = srr, Tissue = '', Genotype = sample_alias,
                  Treatment = '', Replicate = '',
                  paired = paired, spots = spots, avgLength = avgLength)
    #}}}
} else if(yid == 'dn17b') {
    #{{{ hapmap3
    th = ti %>%
        separate("sample_alias", c('pre', 'Genotype'), sep = "_", fill = 'left') %>%
        transmute(SampleID = srr, Tissue = '', Genotype = Genotype,
                  Treatment = '', Replicate = '',
                  paired = paired, spots = spots, avgLength = avgLength)
    #}}}
} else if(yid == 'm282') {
    #{{{ 282set
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = str_sub(SampleName,8),
                  Treatment = '',
                  Replicate = '',
                  paired = paired, spots = spots, avgLength = avgLength) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'gem31') {
    #{{{ german31
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = str_sub(SampleName,4),
                  Treatment = '',
                  Replicate = '',
                  paired = paired, spots = spots, avgLength = avgLength) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'cs845') {
    #{{{
    th = ti %>%
        mutate(SampleName = str_replace_all(SampleName, "^(CRT|PGT)-(\\d+)-(\\d+)_", "\\1_\\2_\\3-")) %>%
        mutate(SampleName = str_replace_all(SampleName, "^(CRT|PGT)-(\\d+)_", "\\1_\\2-")) %>%
        separate(SampleName, c('idx','SampleName'), sep='-', extra='merge') %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = idx,
                  Treatment = SampleName,
                  Replicate = '',
                  paired = paired, spots = spots, avgLength = avgLength) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'me10a') {
#{{{ Li2010
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = 'leaf',
                  Genotype = 'B73',
                  Treatment = SampleName,
                  Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me11a') {
#{{{ Davidson2011
    tismap = c(
"Embryo 25 days after pollination" = 'embryo_25DAP',
"Endosperm 25 days after pollination" = 'endosperm_25DAP',
"Leaves 20-day old seedling - field" = 'leaf_20d_field',
"Leaves 20-day old seedling - growth chamber" = 'leaf_20d_gc',
"Mature silk" = 'silk',
"Ovule" = 'ovule',
"Pollen" = 'pollen',
"Post-emergence cob" = 'cob_post-emerg',
"Postemergence tassel" = 'tassel_post-emerg',
"Pre-emergence cob" = 'cob_pre-emerg',
"Preemergence tassel" = 'tassel_pre-emerg',
"Seed 10 days after pollination" = 'seed_10DAP',
"Seed 5 days after pollination" = 'seed_5DAP',
"Whole anthers" = 'anther')
    th = ti %>%
        separate(Title, c('pre','tis'), sep = ' B73 ') %>%
        separate(tis, c('tis', 'suf'), sep = ' RNA-Seq ') %>%
        mutate(Tissue = tismap[tis]) %>%
        separate(Tissue, c('Tissue','Treatment'), extra='merge', fill='right') %>%
        transmute(SampleID = Run,
                  Tissue = Tissue,
                  Genotype = 'B73',
                  Treatment = Treatment,
                  Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me13a') {
#{{{ Li2013
    th = ti %>% separate("SampleName", c('org', 'ibm', 'parent', 'tis1', 'tis2', 'Genotype'), sep = "_", fill = 'left')
    th %>% count(parent)
    th %>% count(tis1)
    th %>% count(tis2)
    th %>% count(Genotype)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'SAM',
                          Genotype = Genotype,
                          Treatment = '',
                          Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me13b') {
#{{{ Liu2013
    th = ti %>%
        separate("SampleName", c("pre", "Treatment"), sep = "_", fill = "left") %>%
        mutate(Treatment = ifelse(is.na(pre), Treatment, sprintf("E%s", Treatment)))
    th %>% count(paired)
    th %>% count(Treatment)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'leaf',
                          Genotype = 'B73',
                          Treatment = Treatment,
                          Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me13c') {
#{{{ Eitchen2013
    th = ti %>%
        mutate(gt = str_replace(SampleName, "[ _](rep|R) ?[0-9]+", '')) %>%
        mutate(gt = ifelse(gt=='M37W','M37w',gt))
    th %>% count(paired)
    th %>% count(gt)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'unknown',
                          Genotype = gt,
                          Treatment = '',
                          Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me13d') {
#{{{ Waters2013
    th = ti %>% transmute(SampleID = Run,
                          Tissue = 'endosperm',
                          Genotype = SampleName,
                          Treatment = '',
                          Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
    th %>% count(Genotype)
    th %>% count(paired)
#}}}
} else if (yid == 'me13e') {
#{{{ Fu2013
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'kernel',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me14a') {
#{{{ Hirsch2014
th = ti %>% separate(Title, c("org", "Genotype"), sep = ", ") %>%
    separate(Genotype, c("Genotype", "suf"), sep = " RNAseq") %>%
    select(-org, -suf) %>%
    transmute(SampleID = Run,
              Tissue = "seedling",
              Genotype = Genotype,
              Treatment = '',
              Replicate = 1,
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me14b') {
#{{{ Li2014 endosperm
th = ti %>% 
    separate(Title, c("pre", "Treatment"), sep = " B73 ") %>%
    select(-pre) %>%
    separate(Treatment, c("Treatment", "suf"), sep = "DAP; ") %>%
    select(-suf) %>%
    mutate(Treatment = ifelse(Treatment %in% c("0a", "0b"), "0", Treatment)) %>%
    mutate(Treatment = as.integer(Treatment)) %>%
    transmute(SampleID = Run,
              Tissue = "endosperm",
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me14c') {
#{{{ Chettoor gamete 7
th = ti %>%
    mutate(tisrep = str_replace(LibraryName, 'Maize', '')) %>%
    mutate(tis = str_to_lower(str_sub(tisrep, 1, -2)),
           rep = as.integer(str_sub(tisrep, -1, -1))) %>%
    mutate(tis = ifelse(tis == 'embryosac', 'embryo_sac', tis)) %>%
    separate(tis, c("tis","treat"), sep="_", fill='right',extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = tis,
              Genotype = 'B73',
              Treatment = treat,
              Replicate = rep,
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me14d') {
#{{{ johnston ligule
th = ti %>% separate(Title, c("gsm1", "gsm"), sep = ": ") %>%
    separate(gsm, c("tisrep", 'zm', 'rna'), sep = "; ") %>%
    mutate(tisrep = str_replace(tisrep, '[\\.]', '-')) %>%
    separate(tisrep, c('tis', 'rep'), sep = '-') %>%
    transmute(SampleID = Run,
              Tissue = 'ligule',
              Genotype = 'B73',
              Treatment = tis,
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me14e') {
#{{{ Chen seed
th = ti %>%
    mutate(SampleName=str_replace(SampleName, "_B73$", '')) %>%
    mutate(SampleName=str_replace(SampleName, "-DAP", "DAP")) %>%
    mutate(SampleName=str_replace(SampleName, "^embryo_(\\d+)-?DAP$", '\\1DAP_embryo')) %>%
    separate(SampleName, c("stage","tis"), sep="_", fill='right') %>%
    mutate(tis=str_replace(tis, '^whole-','')) %>%
    mutate(tis=str_replace(tis, '^endopserm$', 'endosperm')) %>%
    mutate(stage=str_replace(stage, '^(\\d+)$', '\\1DAP')) %>%
    transmute(SampleID = Run,
              Tissue = tis,
              Genotype = 'B73',
              Treatment = stage,
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me15a') {
#{{{ Leiboff2015
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'SAM',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me15b') {
#{{{ Yu2015
    th = ti %>% separate(Title, c("pre", "Treatment"), sep = " at ") %>%
        transmute(SampleID = Run,
                  Tissue = "leaf",
                  Genotype = 'B73',
                  Treatment = Treatment,
                  Replicate = 1,
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        filter(paired) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'me16b') {
#{{{ Stelpflug2016
th = ti %>%
    mutate(Title = str_replace(Title, 'RNA-seq', 'RNA-Seq')) %>%
    separate(Title, c('pre', 'str1'), sep = ", B73 ", fill='left') %>%
    mutate(str1 = str_replace(str1, '^Zea mays ', '')) %>%
    mutate(str1 = str_replace(str1, ' \\(.*\\)$', '')) %>%
    separate(str1, c('tis_str', 'suf'), sep = " RNA-Seq", fill='right') %>%
    mutate(tis_str = str_replace(tis_str, '_R[1-3]$', '')) %>%
    mutate(tis_str = str_replace(tis_str, ' Rep[1-3]$', ''))
tmap = th %>% select(SampleID=Run, Tissue=tis_str) %>% 
    mutate(nTissue = '') %>% count(Tissue, nTissue)
fo = file.path(dird, '05_read_list/me16b_map_raw.tsv')
write_tsv(tmap, fo)
fo = file.path(dird, '05_read_list/me16b_map.tsv')
tmap = read_tsv(fo)
th = th %>% 
    inner_join(tmap, by = c('tis_str' = 'Tissue')) %>%
    separate(nTissue, c("Tissue",'Treatment'), sep='_', fill='right', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = Tissue,
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me16c') {
#{{{ Walley2016
tismap = c(
"2-4 mm from tip of ear primordium" = 'ear_2-4',
"6-8 mm from tip of ear primordium" = 'ear_6-8',
"Cortex" = 'root_cortex',
"EMBRYOS" = 'embryo',
"embryos_20DAP" = 'embryo_20DAP',
"endosperm_12DAP" = 'endosperm_12DAP',
"endosperm_crown" = 'endosperm_crown',
"EZ" = 'root_ez',
"Germinating Kernels" = 'kernel_germinating',
"GROWTH ZONE" = 'leaf_growth',
"Internode 6-7" = 'internode_6-7',
"Internode 7-8" = 'internode_7-8',
"mature female spikelets" = 'spikelet',
"MATURE LEAF TISSUE (leaf 8)" = 'leaf_mature_8',
"Mature pollen" = 'pollen',
"MZ" = 'root_mz',
"pericarp_aleurone" = 'seed_pericarp',
"PR" = 'root_primary',
"silks" = 'silk',
"SR" = 'root_secondary',
"STOMATAL_DIVISION_ZONE" = 'leaf_stomatal',
"SYMMETRICAL_DIVISION_ZONE" = 'leaf_symmetrical',
"Vegetative Meristem & Surrounding Tissue" = 'meristem')
th = ti %>% separate(Title, c('gsm1', 'gsm'), sep = ': ') %>%
    select(-gsm1) %>%
    separate(gsm, c('tisrep', 'suf1', 'suf2'), sep = '; ') %>%
    select(-suf1, -suf2) %>%
    mutate(tisrep = str_replace(tisrep, "_r([1-3])$", "=\\1")) %>%
    mutate(tisrep = str_replace(tisrep, "_rep([1-3])$", "=\\1")) %>%
    mutate(tisrep = str_replace(tisrep, " ([1-3])$", "=\\1")) %>%
    separate(tisrep, c('Tissue', 'Replicate'), sep = "=") %>%
    mutate(Tissue = tismap[Tissue]) %>%
    separate(Tissue, c("Tissue","Treatment"), sep='_', fill='right') %>%
    transmute(SampleID = Run,
              Tissue = Tissue,
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = Replicate,
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th %>% count(Tissue) %>% print(n=23)
#}}}
} else if (yid == 'me17a') {
#{{{ Lin2017
th1 = ti %>% filter(paired) %>%
    mutate(LibraryName = ifelse(LibraryName=='Mo18W', 'Mo18W-root', LibraryName)) %>%
    separate(LibraryName, c('gt','tissue'), sep = "-") %>%
    mutate(tissue = ifelse(tissue == 'fieldear', 'ear', tissue)) %>%
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = gt,
              Treatment = '',
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th1 %>% count(Tissue)
#
th2 = ti %>% filter(!paired) %>%
    transmute(SampleID = Run,
              Tissue = 'SAM',
              Genotype = LibraryName,
              Treatment = '',
              Replicate = '',
                  paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = rbind(th1,th2)
th %>% count(Tissue)
th %>% count(Genotype)
#}}}
} else if (yid == 'me17c') {
#{{{ Marcon2017
th = ti %>% separate(SampleName, c("gt",'trea','rep'), by='-') %>%
    transmute(SampleID = Run, Tissue = 'root',
              Genotype = gt,
              Treatment = trea,
              Replicate = rep,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me18a') {
#{{{ Kremling2018
th1 = ti %>%
    separate("LibraryName", c('lib1', 'lib2', 'tissue', 'genotype', 'suf'),
             sep = "_", fill = 'left', extra = 'merge')
th1 %>% count(tissue)
tissues = "GRoot GShoot Kern L3Base L3Mid L3Tip LMAD26 LMAD8 LMAN26 LMAN8 LMid"
tissues = strsplit(tissues, " ")[[1]]
idxs = (!th1$tissue %in% tissues)
#
th2 = th1 %>%
    mutate(
           genotype = ifelse(idxs, sprintf("%s_%s", tissue, genotype), genotype),
           tissue = ifelse(idxs, lib2, tissue),
           lib2 = ifelse(idxs, '', lib2))
th2 %>% count(tissue)
#
th3 = th2 %>%
    mutate(tissue = ifelse(tissue %in% c("LMAD26", "LMAD8"), "LMAD", tissue)) %>%
    mutate(tissue = ifelse(tissue %in% c("LMAN26", "LMAN8"), "LMAN", tissue))
th3 %>% count(tissue)
#
tissues = "GRoot GShoot Kern L3Base L3Tip LMAD LMAN"
tissues = strsplit(tissues, " ")[[1]]
th4 = th3 %>% filter(tissue %in% tissues) %>%
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = genotype,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#
th = th4 %>% mutate(Treatment = '', Replicate = '') %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, spots, avgLength)
th %>% count(Tissue) %>% print(n=100)
#}}}
} else if (yid == 'me18b') {
#{{{ Baldauf2018
th = ti %>%
    mutate(SampleName = str_replace(SampleName, '-', '_0_')) %>%
    separate(SampleName, c("gt", "stage", "rep"), sep = "_") %>%
    transmute(SampleID = Run,
              Tissue = 'root',
              Genotype = gt,
              Treatment = stage,
              Replicate = rep,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me18c') {
#{{{ Wang2018
fd = file.path(dird, '03_sra_list','me18c.txt')
td = read_tsv(fd, col_names=c('yid','sid2'))
th = ti %>%
    transmute(SampleID = Run,
              Tissue = 'seedling_v2',
              Genotype = LibraryName,
              Treatment = NA,
              Replicate = NA,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    filter(Genotype %in% td$yid) %>%
    arrange(SampleID)
#}}}
} else if (yid == 'me18d') {
#{{{ Schaefer2018
th = ti %>%
    separate(SampleName, c("gt", "suf1", "suf2"), sep="-", fill='right', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = 'root',
              Genotype = gt,
              Treatment = '',
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'me18e') {
#{{{ Huang2018
    fh = sprintf("%s/05_read_list/%s_raw.tsv", dird, yid)
    tiss = c("root",'leaf','SAM','seed')
    th2 = read_tsv(fh) %>%
        mutate(Tissue=ifelse(Tissue %in% tiss, Tissue, 'seed')) %>%
        select(SampleID, Tissue, Genotype)
th = ti %>%
    select(SampleID=Run, Treatment=Experiment,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    inner_join(th2, by = 'SampleID') %>%
    mutate(Replicate = '') %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired,spots,avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'me19a') {
#{{{
th = ti %>%
    transmute(SampleID = Run,
              Tissue = "seedling",
              Genotype = SampleName,
              Treatment = '',
              Replicate = 1,
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'me99a') {
#{{{ Kaeppler2018
th = ti %>%
    mutate(gt0 = str_replace(Title2, "^.*Zea mays ?", "")) %>%
    mutate(gt0 = str_replace(gt0, " ?transcriptome$", '')) %>%
    mutate(gt0 = str_replace(gt0, " ?gene expression profiling.*$", '')) %>%
    mutate(gt0 = str_replace(gt0, ' ', '')) %>%
    mutate(gt0 = str_replace(gt0, "_([0-9]{4})$", "~\\1")) %>%
    separate(gt0, c("gt1",'tis'), sep='_', fill='right',extra='merge') %>%
    mutate(gt1 = str_replace(gt1, "~([0-9]{4})$", "_\\1")) %>%
    replace_na(list(tis='S')) %>%
    mutate(tis=str_replace(tis,'_T2', '')) %>%
    mutate(tis=str_replace(tis,'^L\\?\\?$','L')) %>%
    transmute(SampleID = Run, Tissue = tis, Genotype = gt1,
              Treatment = '', Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
idx = which(th$Tissue=='MoG_115')
th$Genotype[idx] = sprintf("%s-MoG-115", th$Genotype[idx])
th$Tissue[idx] = 'S'
tismap = c('S'='seedling','I'='internode','R'='root','L'='leaf','E'='endosperm')
th = th %>%
    mutate(Tissue=tismap[Tissue]) %>%
    mutate(inbred = !str_detect(Genotype, '[xX]'))
th %>% count(Treatment) %>% print(n=20)
th %>% filter(inbred) %>% distinct(Tissue)
th %>% filter(inbred) %>% distinct(Genotype) %>% pull(Genotype)
th %>% filter(!inbred) %>% count(Tissue)
th %>% filter(!inbred) %>% distinct(Genotype) %>% pull(Genotype)
#}}}
} else if (yid == 'mem01') {
#{{{ Bolduc2012
    th = ti %>%
        separate(LibraryName, c("gsm1", "gsm"), sep = ": ") %>%
        separate(gsm, c("gentisrep", 'zm', 'rna'), sep = "; ") %>%
        separate(gentisrep, c("gentis", 'rep'), sep = " #") %>%
        mutate(gentis = str_replace(gentis, " leaf (homo|het)", "_\\1 leaf")) %>%
        separate(gentis, c("gen", "tis"), sep = " ") %>%
        mutate(tis = str_replace(tis, 's$', '')) %>%
        transmute(SampleID = Run,
                  Tissue = tis,
                  Genotype = 'B73',
                  Treatment = gen,
                  Replicate = rep,
                  paired = paired, spots = spots, avgLength=avgLength) %>%
        arrange(SampleID)
#}}}
} else if (yid == 'mem02') {
#{{{
th = ti %>% separate(LibraryName, c('tre','rep'), sep='_') %>%
    transmute(SampleID = Run,
              Tissue = 'ear tip',
              Genotype = 'B73',
              Treatment = tre,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem03') {
#{{{
th = ti %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ') %>%
    mutate(tre = str_replace(tre, " Rep[0-9]$", "")) %>%
    transmute(SampleID = Run,
              Tissue = 'endosperm',
              Genotype = 'B73',
              Treatment = tre,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem04') {
#{{{
th = ti %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    mutate(tre = ifelse(tre == 'wild type', 'WT', 'mads47')) %>%
    transmute(SampleID = Run,
              Tissue = 'kernel',
              Genotype = 'B73',
              Treatment = tre,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem05') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('gt','rep','tis'), sep='_', extra='merge') %>%
    mutate(gt = ifelse(gt %in% c('wt','B73'), 'WT', gt)) %>%
    transmute(SampleID = Run,
              Tissue = tis,
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem06') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('tis','gt','rep'), sep='_', extra='merge') %>%
    mutate(gt = ifelse(gt %in% c('wildtype'), 'WT', 're2')) %>%
    transmute(SampleID = Run,
              Tissue = 'tassel',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem07') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('gt','suf','rep'), sep=' ', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = 'endosperm',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem08') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    mutate(gt = ifelse(tre == 'wild type', 'WT', 'o2')) %>%
    transmute(SampleID = Run,
              Tissue = 'endosperm',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem09') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('gt','rep'), sep='_', extra='merge') %>%
    mutate(gt = ifelse(gt == 'zmbzip22', 'bzip22', gt)) %>%
    transmute(SampleID = Run,
              Tissue = 'kernel',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem10') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('gt','rep'), sep='-', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = 'kernel',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem11') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('gt','rep'), sep='_', extra='merge') %>%
    mutate(gt = ifelse(gt == 'fea4', gt, 'WT')) %>%
    transmute(SampleID = Run,
              Tissue = 'ear',
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem12') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = 'leaf primordia',
              Genotype = 'B73',
              Treatment = tre,
              Replicate = '',
              paired = paired, spots = spots, avgLength=avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else if (yid == 'mem13') {
#{{{
th = ti %>% filter(LibraryStrategy == 'RNA-Seq') %>%
    separate(Title, c('str1','str2'), sep=': ') %>%
    separate(str2, c('tre','str2b','str2c'), sep='; ', extra='merge') %>%
    separate(tre, c('tre','str3'), sep=', ', extra='merge') %>%
    separate(tre, c('gt','tis'), sep=' ', extra='merge') %>%
    mutate(tis = ifelse(tis == 'Aleurone', 'aleurone', 'endosperm')) %>%
    mutate(gt = ifelse(gt == 'B73', 'WT', gt)) %>%
    transmute(SampleID = Run,
              Tissue = tis,
              Genotype = 'B73',
              Treatment = gt,
              Replicate = '',
              paired = paired, spots = spots, avgLength = avgLength) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
#}}}
} else {
    cat("unknown study: ", yid, "\n")
}
th %>% arrange(SampleID) %>% sra_fill_replicate(th)
#}}}
}


