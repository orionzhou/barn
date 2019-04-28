source("functions.R")
t_cfg %>% select(-yid2,-stranded,-interleaved) %>% print(n=40)

yid = 'sp052'
opt = t_cfg %>% filter(yid==!!yid) %>% pull(source)
interleaved = t_cfg %>% filter(yid==!!yid) %>% pull(interleaved)
fi = file.path(dird, '03.raw.xlsx')
ti = read_xlsx(fi, sheet=yid)
to = complete_sample_list(ti) %>%
    mutate(data = map2(directory,file_prefix, locate_fastq,
                       opt=!!opt, interleaved=!!interleaved)) %>%
    unnest() %>% select(-directory,-file_prefix) %>%
    group_by(Tissue,Genotype,Treatment) %>%
    mutate(Replicate = 1:n()) %>%
    ungroup()
to %>% count(Genotype, Tissue, Treatment) %>% print(n=50)
#
fo = sprintf("%s/10_raw_list/%s.tsv", dird, yid)
write_tsv(to, fo)



