## Run from BISCUT dev directory

library(devtools)
load_all() # load dev version of package
library(data.table)

# Get pan-TCGA peaks from Shih et al., table S3A
original = as.data.table(readxl::read_excel('for_compare/41586_2023_6266_MOESM5_ESM.xlsx', 
                                            sheet = 1, skip = 17))
original = original[, peak_id := paste(arm, direction, telcent, code, iter, sep = "_")]


# Compare to peaks running on altered version of BISCUT-py3
# Downloaded from the link given in the BISCUT README: 
# https://drive.google.com/drive/folders/1-VZ_A0uodEs04Jg-Gphkl3HU5VdEyICW?usp=share_link
segment_file = 'for_compare/PANCAN_ISAR.seg.txt.gz'

# You may need to adjust cores depending on your system
bp = calculate_breakpoints(segment_file = segment_file, cores = 8)
results = do_biscut(breakpoints = bp, cores = 8, seed = 999)

original_peaks = unique(original[, .(peak_id, code, iter, direction, telcent, negpos, combined_sig)])
new_peaks = results$peaks[, .(peak_id, code, iter, direction, telcent, negpos, combined_sig, search_lowlim, search_highlim,
                              Peak.Start.1, Peak.End.1, log10_ksby)]


all_peak_id = unique(c(original_peaks$peak_id, new_peaks$peak_id))
peak_compare = data.table(peak_id = all_peak_id,
                          in_original = all_peak_id %in% original_peaks$peak_id,
                          in_new = all_peak_id %in% new_peaks$peak_id)


# 92.7% of the peaks in the original results match peaks found in the new run.
peak_compare[in_original == T, mean(in_new)]
# [1] 0.9274611

table(peak_compare[, .(in_original, in_new)])
#             in_new
# in_original FALSE TRUE
# FALSE       0   41
# TRUE        14  179

# Considering only 1st iteration peaks, 97.9% of original peaks are reproduced in the new results.
peak_compare[, iter := sub('.*_', '', peak_id)]
first_iter_compare = peak_compare[iter == 1]
first_iter_compare[in_original == T, mean(in_new)]
# [1] 0.9789474


table(first_iter_compare[, .(in_original, in_new)])
#             in_new
# in_original FALSE TRUE
# FALSE       0   7
# TRUE        2   93

original_gene_results = original[, .(gene = Gene, negpos, ampdel = direction, telcent)]
orig_single_hit = original_gene_results[, .N, by = names(original_gene_results)][N == 1, gene]
orig_for_compare = original_gene_results[gene %in% orig_single_hit]

curr_gene_results = results$genes_by_peak[, .(peak_id, gene = Gene)]
curr_gene_results[results$peaks, let(negpos = negpos, ampdel = direction, telcent = telcent), on = 'peak_id']
curr_single_hit = curr_gene_results[ ,.N, by = names(curr_gene_results)][N == 1, gene]
curr_for_compare= curr_gene_results[gene %in% curr_single_hit]

# Compare the 8,531 gene-ampdel-telcent peaks that are shared between the runs.
gene_compare = merge.data.table(orig_for_compare, curr_for_compare, suffixes = c('.original', '.new'), 
                                by = c('gene', 'ampdel', 'telcent'), all = F)
stopifnot(gene_compare[, .N] == 8531)

## We get 100% agreement!
gene_compare[, mean(negpos.original == negpos.new)]
# [1] 1



