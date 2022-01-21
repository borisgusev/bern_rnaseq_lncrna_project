library('sleuth')
library('tidyselect')
library('dplyr')
library('ggplot2')


# load data in format accepted by sleuth
sample_num <- seq(0, 11)
res_dirs <- file.path('..', 'samples', sample_num)
metadata <- read.table(file.path('..', 'samples', 'metadata.txt'), header = T, stringsAsFactors = T)
metadata$condition <- relevel(metadata$condition, 'parental')
metadata <- dplyr::mutate(metadata, path = res_dirs)


# load the transcript -> gene mapping
t2g <- read.csv('genename_gid_tid.tsv', header = F, stringsAsFactors = F, sep = '\t')
colnames(t2g) <- c('gene_name', 'gene_id', 'target_id')

#--------------------------
# Transcript Level Analysis
#--------------------------

# make and manipulate sleuth object as per the tutorial 
# added function to have log2 instead of the default natural log
so <- sleuth_prep(metadata, extra_bootstrap_summary = T, read_bootstrap_tpm = T, transformation_function = function(x) log2(x + 0.5))
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')


# assemble the transcript level expression table
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
qvals <- dplyr::select(sleuth_table, target_id, qval)

transcript_level_df <- tibble::as_tibble(sleuth_to_matrix(so, which_df = 'obs_norm', which_units = 'tpm'), rownames = 'target_id')
transcript_level_df <- dplyr::inner_join(transcript_level_df, qvals)
transcript_level_df <- dplyr::rename(transcript_level_df, transcript_id = target_id)
readr::write_csv(transcript_level_df, 'transcript_level_expression_data.csv')


# plotting

#pca
plot_pca(so, color_by = 'condition') + 
  ggtitle('PCA of clonal expression levels')
ggsave('pca.png', path = 'plots/')

# volcano plots for each condition using wald test

so <- sleuth_wt(so, 'conditionholoclone', 'full')
plot_volcano(so, 'conditionholoclone') +
  ggtitle('Holoclone Differential Expression Volcano Plot')
ggsave('volcano_holoclone.png', path = 'plots/')


so <- sleuth_wt(so, 'conditionmeroclone', 'full')
plot_volcano(so, 'conditionmeroclone') + 
  ggtitle('Meroclone Differential Expression Volcano Plot')
ggsave('volcano_meroclone.png', path = 'plots/')

so <- sleuth_wt(so, 'conditionparaclone', 'full')
plot_volcano(so, 'conditionparaclone') + 
  ggtitle('Paraclone Differential Expression Volcano Plot')
ggsave('volcano_paraclone.png', path = 'plots/')

#sox2
plot_bootstrap(so, "ENST00000325404.3", units = "tpm", color_by = "condition") +
  ggtitle('SOX2 Transcript Expression')
ggsave('expression_sox2.png', path = 'plots/')
#cdh1
plot_bootstrap(so, "ENST00000261769.10", units = "tpm", color_by = "condition") +
  ggtitle('CDH1 Transcript Expression')
ggsave('expression_cdh1.png', path = 'plots/')
# 
# plot_bootstrap(so, "ENST00000257555.11", units = "tpm", color_by = "condition")
# 





#---------------------#
# Gene Level Analysis #
#---------------------#

so_gene <- sleuth_prep(metadata, 
                       extra_bootstrap_summary = T, 
                       read_bootstrap_tpm = T, 
                       target_mapping = t2g, 
                       gene_mode = T, 
                       aggregation_column = 'gene_id', 
                       transformation_function = function(x) log2(x + 0.5))

so_gene <- sleuth_fit(so_gene, ~condition, 'full')
so_gene <- sleuth_fit(so_gene, ~1, 'reduced')
so_gene <- sleuth_lrt(so_gene, 'reduced', 'full')

# assemble the transcript level expression table
sleuth_table_gene <- sleuth_results(so_gene, 'reduced:full', 'lrt')
qvals_gene <- dplyr::select(sleuth_table_gene, target_id, gene_name, qval)

gene_level_df <- tibble::as_tibble(sleuth_to_matrix(so_gene, which_df = 'obs_norm', which_units = 'tpm'), rownames = 'target_id')
gene_level_df <- dplyr::inner_join(gene_level_df, qvals_gene)
gene_level_df <- dplyr::rename(gene_level_df, gene_id = target_id)
readr::write_csv(gene_level_df, 'gene_level_data.csv')



sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
qvals <- dplyr::select(sleuth_table, target_id, qval)

transcript_level_df <- tibble::as_tibble(sleuth_to_matrix(so, which_df = 'obs_norm', which_units = 'tpm'), rownames = 'target_id')
transcript_level_df <- dplyr::inner_join(transcript_level_df, qvals)
transcript_level_df <- dplyr::rename(transcript_level_df, transcript_id = target_id)
readr::write_csv(transcript_level_df, 'transcript_level_expression_data.csv')



# # plot_bootstrap(so_gene, "MSTRG.21157", units = "tpm", color_by = "condition")






