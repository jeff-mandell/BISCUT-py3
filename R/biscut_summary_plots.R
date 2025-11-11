#' Produce peak plots
#' 
#' Visualize significant peaks
#' 
#' @param biscut_results The data (list) that is output from do_biscut()
#' @param plot_dir Path for a new directory to save plot PDFs (or NULL to not make any).
#' @return A list containing all the summary plots.
#' @export
biscut_summary_plots = function(biscut_results, plot_dir = NULL) {
  if(! is.null(plot_dir)) {
    if(! is.character(plot_dir) || length(plot_dir) != 1 || nchar(plot_dir) == 0) {
      stop('plot_dir should be scalar character (path for plot output directory), or left NULL.')
    }
    if(! dir.exists(dirname(plot_dir))) {
      stop('Parent directory for input plot_dir (', dirname(plot_dir), ' does not exist.')
    }
    if(dir.exists(plot_dir)) {
      stop('Specified plot_dir directory (', plot_dir, ') already exists.')
    }
    if(! dir.create(plot_dir)) {
      stop('Could not create plot_dir')
    }
    plot_dir = normalizePath(plot_dir)
  }
  if(! is.list(biscut_results) || ! 'peaks' %in% names(biscut_results)) {
    stop('biscut_results should be type list (the output from do_biscut()).')
  }
  
  abslocs = copy(biscut_results$info$abslocs)
  abslocs_dt = as.data.table(abslocs) # for simpler merging into peak data
  if(! 'offset' %in% names(abslocs_dt)) {
    abslocs_dt$size = as.numeric(abslocs_dt$size) # avoid integer overflow (!)
    abslocs_dt[, offset := c(0, cumsum(abslocs_dt$size[1:(.N - 1)]))]
  }
  
  
  for_plots = copy(biscut_results$peaks)
  for_plots[direction == 'del', log10_ksby := -1 * log10_ksby]
  for_plots[, combined_sig := log10_ksby * ks_stat]
  
  darkred = '#a50f15'
  lightred = '#fcae91'
  lightblue = '#6baed6'
  darkblue = '#08519c'
  
  # Merge in offset (effectively converting chromosome coordinates to a genome-wide coordinate system)
  for_plots[abslocs_dt, offset := offset, on = c(Chr = 'chromosome_info')]
  for_plots[, xmin := Peak.Start + offset]
  for_plots[, xmax := Peak.End + offset]
  for_plots[, ymin := pmin(0, combined_sig)]
  for_plots[, ymax := pmax(0, combined_sig)]
  
  write_plot_data = function(dt, name, neg_color, pos_color) {
    if(dt[, .N] == 0) {
      return(data.table())
    }
    dt[pmin(ymin, ymax) < 0, color := neg_color]
    dt[pmax(ymin, ymax) > 0, color := pos_color]
    output = dt[, .(peak_id, direction, type_of_selection, xmin, xmax, ymin, ymax, color)]
    
    if(! is.null(plot_dir)) {
      output_file = paste0(plot_dir, '/plot_data_', name, '.txt')
      fwrite(output, file = output_file, sep = "\t")
    }
    return(output)
  }
  
  amp_data = write_plot_data(for_plots[direction == 'amp'], name = 'amp', neg_color = lightblue, pos_color = darkred)
  del_data = write_plot_data(for_plots[direction == 'del'], name = 'del', neg_color = lightred, pos_color = darkblue)
  pos_data = write_plot_data(for_plots[negpos == 'p'], name = 'pos', neg_color = darkblue, pos_color = darkred)
  neg_data = write_plot_data(for_plots[negpos == 'n'], name = 'neg', neg_color = lightred, pos_color = lightblue)
  
  xlab1 = "Significance\nnegative selection <----------------------------------------------> positive selection"
  amp_plot = plot_BISCUT_results(amp_data, title = 'Amplifications', xlab = xlab1, chromosome_coordinates = abslocs)
  del_plot = plot_BISCUT_results(del_data, title = 'Deletions', xlab = xlab1, chromosome_coordinates  = abslocs)
  
  xlab2 = "Significance\nessential-like <----------------------------------------------> toxic-like"
  neg_plot = plot_BISCUT_results(neg_data, title = 'Negative Selection', xlab = xlab2, chromosome_coordinates = abslocs)
  
  xlab3 = "Significance\ntumor-suppressor-like <----------------------------------------------> oncogene-like"
  pos_plot = plot_BISCUT_results(pos_data, title = 'Positive Selection', xlab = xlab3, chromosome_coordinates = abslocs)
  
  if(! is.null(plot_dir)) {
    if(amp_data[, .N] > 0) {
      cowplot::save_plot(filename = paste0(plot_dir, '/BISCUT_peaks_', 'amplifications.pdf'), amp_plot,
                         base_width = 8, base_height = 6)
    }
    if(del_data[, .N] > 0) {
      cowplot::save_plot(filename = paste0(plot_dir, '/BISCUT_peaks_', 'deletions.pdf'), del_plot,
                         base_width = 8, base_height = 6)
    }
    if(pos_data[, .N] > 0) {
      cowplot::save_plot(filename = paste0(plot_dir, '/BISCUT_peaks_', 'positive_selection.pdf'), pos_plot,
                         base_width = 8, base_height = 6)
    }
    if(neg_data[, .N] > 0) {
      cowplot::save_plot(filename = paste0(plot_dir, '/BISCUT_peaks_', 'negative_selection.pdf'), neg_plot,
                         base_width = 8, base_height = 6)
    }
  }
  return(list(amplifications = amp_plot, deletions = del_plot,
              positive_selection = pos_plot, negative_selection = neg_plot))
}


