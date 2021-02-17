source("functions.R")

genome = 'at'
genome = 'zm'
genome = 'os'

t_cfg = read_projects(genome)
t_cfg %>% count(libtype,source)
t_cfg %>% print(n=50)
rcfg = t_cfg %>% filter(source=='sra')
lcfg = t_cfg %>% filter(source=='local')

yid = 'mo20b'
acc = rcfg %>% filter(yid == !!yid) %>% pull(accession)
ti = get_sra_meta(acc, yid)
#
fo = glue("{dird}/tmp.tsv")
write_tsv(ti, fo)

#{{{ SRA - work on single one
yid = 'cg18a'
acc = rcfg %>% filter(yid == !!yid) %>% pull(accession)
if(yid == 'rn99a') acc = read_tsv(file.path(dird, 'rn99a_bpid.txt'),col_names=F)$X1
ti = get_sra_meta(acc, yid)

fo = glue("{dirw}/tmp.tsv")
write_tsv(ti, fo)

#tx = read_xlsx(fx)
#tx = complete_sample_list(tx, fillrep=T)
#fo = sprintf("%s/09_sra_list/%s.tsv", dird, yid)
#write_tsv(tx, fo, na = '')
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
yid = 'rn20h'
lid = lcfg %>% filter(yid==!!yid) %>% pull(accession)
fi = sprintf("%s/05_msi/%s.tsv", dird, lid)
ti = read_tsv(fi)
to = complete_sample_list(ti, fillrep=F) %>%
    mutate(data = map2(directory,file_prefix, locate_fastq,
                       fmt='umgc', interleaved=F)) %>%
    unnest(data) %>% select(-directory,-file_prefix)
to %>% count(MergeID, Genotype, Tissue, Treatment) %>% print(n=50)

fo = sprintf("%s/06_local_list/%s.tsv", dird, yid)
write_tsv(to, fo, na='')
#}}}




