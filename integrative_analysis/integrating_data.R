library('tidyverse')
library('dplyr')
library('ggplot2')

novel_transcripts <- readr::read_table('cpat/merged.transcripts.novel.gtf', col_names = F) %>%
  select(gene_id = X10, transcript_id = X12) 

cpat_output <- readr::read_tsv('cpat/cpat_output.ORF_prob.tsv') %>%
  filter(Coding_prob < 0.364) %>%
  select(ID) %>%
  filter(endsWith(ID, 'ORF_1')) %>%
  tidyr::separate(ID, c('gene_id', NA), sep = '::') %>%
  mutate(gene_id = paste0('"', gene_id, '";'), isNonCoding = T) # lazy hack to match the gtf inputs that have "ID";

cage_transcripts <- readr::read_table('end_analysis/merged.transcripts.novel.gtf.withcage', col_names = F) %>%
  select(gene_id = X10, transcript_id = X12) %>%
  mutate(hasCage = T)

polya_transcripts <- readr::read_table('end_analysis/merged.transcripts.novel.gtf.withpolya', col_names = F) %>%
  select(gene_id = X10, transcript_id = X12) %>%
  mutate(hasPolya = T)

intergenic_transcripts <- readr::read_table('intergenic/merged.transcripts.novel.intergenic.gtf', col_names = F) %>%
  select(gene_id = X10, transcript_id = X12) %>%
  mutate(isIntergenic = T)

multi_exon_transcripts <- readr::read_table('single_exon/merged.transcripts.exoncount.novel.gtf', col_names =F) %>%
  select(gene_id = X10, transcript_id = X12, exon_count = X13) %>%
  mutate(isMultiExon = (exon_count != 1)) %>%
  select(gene_id, transcript_id, isMultiExon)

integrated <- left_join(novel_transcripts, cage_transcripts) %>%
  left_join(polya_transcripts) %>%
  left_join(cpat_output) %>%
  left_join(intergenic_transcripts) %>%
  left_join(multi_exon_transcripts) %>%
  replace(is.na(.), F)


# best_transcripts <- integrated %>% mutate(total = rowSums(c(hasCage, hasPolya, isIntergenic, isNonCoding, isMultiExon)))

integrated$numTrue <- rowSums(integrated[3:length(integrated)])

integrated_summary <- as.data.frame(t(summarise_all(select(integrated, 4:length(integrated) - 1), mean))) %>%
  rownames_to_column() %>%
  rename(test = 1, True = 2)  %>%
  mutate(False = 1 - True) %>%
  tidyr::pivot_longer(!test) %>%
  rename(Result = name)

ggplot(integrated_summary, aes(x = test, y = value, fill = Result)) +
  geom_bar(stat = 'identity', position = 'fill') +
  xlab('') + ylab('Proportion') +
  ggtitle('Integrative Analysis') 
ggsave('int_an.png', path = 'plots/')


lrt_qvals <- readr::read_csv('../expression_quant/differential_expression/transcript_level_expression_data.csv') %>%
  select(transcript_id, qval) %>%
  mutate(transcript_id = paste0('"', transcript_id, '";')) # lazy hack to match the gtf inputs that have "ID";

integrated <- integrated %>%
  left_join(lrt_qvals) 

readr::write_csv(integrated, 'integrated_analysis.csv')


top_candidates <- integrated %>%
  arrange(desc(numTrue), qval) %>%
  head(40)
readr::write_csv(top_candidates, 'top_novel_candidates.csv')

  



  