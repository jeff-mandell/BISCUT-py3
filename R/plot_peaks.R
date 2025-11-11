#' Produce peak plots
#' 
#' Visualize significant peaks
#' 
#' @param biscut_results The data (list) that is output from do_biscut()
#' @param plot_dir Path for a new directory to save plot PDFs (or NULL to not make any).
#' @return A list containing all the peak plot objects.
#' @export
plot_peaks = function(biscut_results, plot_dir = NULL) {
  # Validate plot directory
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

  message("Making peak plots...")
  peak_info = copy(biscut_results$peaks)
  peak_plot_data = biscut_results$peak_plot_data[peak_info$peak_id] # list, one element per peak_id
  ci = biscut_results$info$ci
  all_plots = list()
  for(i in 1:length(peak_plot_data)) {
    curr_info = peak_info[names(peak_plot_data)[i], on = 'peak_id']
    curr_info[, selection := ifelse(negpos == 'n', 'negative', 'positive')]
    curr_info[, ks_p := 10^(-1 * log10_ks_p)]  # log10_ks_p is actually -1 * log10_ks_p
    title = curr_info[, paste0(paste(arm, direction, telcent, selection, code, ci,'iteration', iter, sep=' '),
                               '\nks p-value = ', signif(ks_p, 3))]
    peak_pos = curr_info$Peak.Position.1
    curr_data = peak_plot_data[[i]]
    curr_data$rownamez = 1:nrow(curr_data)
    p1 = ggplot(curr_data, aes(x = percent, y = rownamez)) + geom_point(size = .75) + xlab('pSCNA Length') + ylab('Ranked Tumors') + ggtitle(title) +
      geom_point(aes(x = x, y = rownamez), shape = '.', size = 1) + 
      geom_vline(xintercept = peak_pos) + theme_classic()
    lowlim = curr_info$search_lowlim
    highlim = curr_info$search_highlim
    left_boundary = curr_info$Peak.Start.1
    right_boundary = curr_info$Peak.End.1
    p2 = ggplot(data= curr_data, aes(x = percent, y = dis)) + geom_point(size = .75) +
      xlab('pSCNA Length') + ylab('Distance from background') + 
      scale_x_continuous(limits = c(lowlim, highlim)) + 
      geom_vline(xintercept = peak_pos) +
      geom_vline(xintercept = left_boundary, lty = 'dashed') + 
      geom_vline(xintercept = right_boundary, lty = 'dashed') + theme_classic()
    peak_plots = cowplot::plot_grid(p1, p2, nrow = 2)
    if(! is.null(plot_dir)) {
      plot_file = paste0(plot_dir, '/', 
                         curr_info[, .(paste0(arm, '_', direction, '_', telcent, '_',
                                              code, '_', ci, '_iter', iter, '.pdf'))])
      cowplot::save_plot(plot_file, peak_plots)
    }
    all_plots[[i]] = peak_plots
  }
  names(all_plots) = names(peak_plot_data)
  return(all_plots)
}
