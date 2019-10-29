source("functions.R")
dirw = file.path(dird, '04_prep')

#{{{ rn18i - ch001
yid = 'rn18i'
gmap = c('B'='B73','P'='PH207','W'='W22')
tmap = c('R'='root','SC'='shoot','I'='internode','L10'='leaf10','L'='leaf30',
    'T'='tassel','IE'='ear','A'='anther','En'='endosperm','Em'='embryo')
dirs = c('/home/hirschc1/data_release/umgc/hiseq/180509_D00635_0371_BCC5NGANXX/Hirsch2_Project_001', '/home/hirschc1/data_release/umgc/hiseq/180522_D00635_0375_BCCJL5ANXX/Hirsch2_Project_001')
tf = tibble(directory = dirs) %>%
    mutate(fn = map(directory, list.files)) %>% unnest() %>%
    filter(str_detect(fn, "\\.fastq\\.gz$")) %>%
    mutate(file_prefix = str_replace(fn, '_R1_001.fastq.gz', '')) %>%
    separate(fn, c('pre','suf'), sep="_", extra='merge') %>%
    separate(pre, c('Genotype','Tissue','Replicate'), sep='-') %>%
    mutate(Genotype = gmap[Genotype], Tissue=tmap[Tissue]) %>%
    mutate(Replicate = as.integer(str_replace(Replicate,'R',''))) %>%
    arrange(Genotype, Tissue, Replicate) %>%
    mutate(SampleID = sprintf("s%03d", 1:n())) %>%
    select(SampleID,Genotype,Tissue,Replicate,directory,file_prefix)

fo = sprintf("%s/%s.tsv", dirw, yid)
write_tsv(tf, fo)
#}}}
