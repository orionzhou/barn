source("functions.R")
dird = '~/projects/rnaseq/data'
dirw = file.path(dird, '11_qc', 'rn20a')

diri = '/home/springer/zhoux379/projects/barn/data/cache/rn20a/umi_cnt'
ti = crossing(x=1:12, y=1:8) %>%
    mutate(yl = LETTERS[y]) %>%
    mutate(fn = str_c(yl,x,sep='')) %>%
    rbind(tibble(x=13, y=9, yl='I',fn='undetermined')) %>%
    mutate(fi = sprintf("%s/%s.fastq.tsv", diri, fn)) %>%
    mutate(res = map(fi, read_tsv, col_names=c('umi','cnt'))) %>%
    select(x, y, yl, fn, res) %>%
    unnest()

#{{{ nseq
tp = ti %>% group_by(x,y,fn) %>%
    summarise(cnt = sum(cnt)) %>% ungroup() %>%
    mutate(lab = number(cnt))
mid = (min(tp$cnt) + max(tp$cnt)) / 2
tp = tp %>% mutate(color = ifelse(cnt < mid, 'white','black'))
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=cnt), color='black') +
    geom_text(aes(x,y, label=lab, color=cnt<mid), hjust=1, size=2.5, nudge_x=.4) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle('# reads') +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b01.nseq.pdf')
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ UMI
tp = ti %>% group_by(x,y,fn) %>%
    summarise(cnt = n()) %>% ungroup() %>%
    mutate(lab = number(cnt))
mid = (min(tp$cnt) + max(tp$cnt)) / 2
tp = tp %>% mutate(color = ifelse(cnt < mid, 'white','black'))
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=cnt), color='black') +
    geom_text(aes(x,y, label=lab, color=cnt<mid), hjust=1, size=2.5, nudge_x=.4) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle('# unique UMIs') +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b02.umi.pdf')
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ umi - distri
tp = ti %>% filter(x<13, cnt<=2000) %>%
    mutate(cnti = cut_interval(cnt, n=10)) %>%
    count(x,yl,cnti)
p = ggplot(tp) +
    geom_bar(aes(cnti, n),stat='identity') +
    scale_x_discrete(name='# UMI occurences', expand=expand_scale(mult=c(.05,.05))) +
    scale_y_continuous(expand=expand_scale(mult=c(.05,.05))) +
    facet_grid(yl ~ x) +
    otheme(xtitle=T, xtext=T, ytext=T, legend.pos='none') +
    theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1)) +
    #theme(panel.border = element_blank()) +
    ggtitle('UMI distribution') +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b05.umi.distr.pdf')
ggsave(fo, p, width=12, height=8)
#}}}

fi1 = file.path(diri, '../dge/output.dge.reads.txt')
fi2 = file.path(diri, '../dge/output.dge.umis.txt')
ti1 = read_tsv(fi1)
ti2 = read_tsv(fi2)
#
tu1 = ti1 %>% rename(gid=1) %>% gather(sid, n_read, -gid)
tu2 = ti2 %>% rename(gid=1) %>% gather(sid, n_umi, -gid)
tu = tu1 %>% inner_join(tu2, by=c('gid','sid'))

l2n = 1:8; names(l2n) = LETTERS[l2n]
min_nr = 2; min_nu = 2
min_nr = 10; min_nu = 10
tp = tu %>% group_by(sid) %>%
    summarise(nr = sum(n_read >= min_nr), nu = sum(n_umi >= min_nu)) %>%
    ungroup() %>%
    mutate(yl = str_sub(sid, 1, 1), x = str_sub(sid, 2)) %>%
    mutate(y = as.integer(l2n[yl]), x = as.integer(x))


#{{{ ngene
mid = (min(tp$nr) + max(tp$nr)) / 2
tit=sprintf('# genes with >= %d uniquely mapped reads', min_nr)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nr), color='black') +
    geom_text(aes(x,y, label=number(nr), color=nr<mid), hjust=1, size=3, nudge_x=.2) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b11.nread.pdf')
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ numi
mid = (min(tp$nu) + max(tp$nu)) / 2
tit=sprintf('# genes with >= %d UMIs', min_nu)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nu), color='black') +
    geom_text(aes(x,y, label=number(nu), color=nu<mid), hjust=1, size=3, nudge_x=.2) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b12.numi.pdf')
ggsave(fo, p, width=8, height=5)
#}}}
