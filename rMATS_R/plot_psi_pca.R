library('ggplot2')

## Rscript plot_psi_pca.R pca_psi.tsv pca_psi.png
args <- base::commandArgs(trailingOnly=TRUE)
in_file_name <- args[1]
out_file_name <- args[2]

create_ggplot_theme <- function() {
  return(
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(color='black'),
                   legend.title=ggplot2::element_blank(),
                   legend.text=ggplot2::element_text(size=4))
  )
}

calculate_pca <- function(data) {
  ## data: columns are samples, rows are events
  ## convert so that rows are samples
  transposed <- base::t(data)
  variances <- base::apply(transposed, 2, stats::var)
  transposed <- transposed[, variances != 0]
  transposed <- base::log10(transposed + 1)
  pca_results <- stats::prcomp(transposed, center=TRUE, scale.=TRUE)
  variance_by_component <- base::apply(pca_results$x, 2, stats::var)
  total_variance <- base::sum(variance_by_component)
  percent_variance <- base::round(
      100 * (variance_by_component / total_variance), digits=1)
  return(list(pca=pca_results, percent_variance=percent_variance))
}

plot_pc_1_2 <- function(plot_df, pca_variance, out_path) {
  point_size <- 1
  point_alpha <- 0.8
  width <- 5
  height<- 4

  title <- 'PCA by PSI value'
  plot <- ggplot2::ggplot(data=plot_df,
                          ggplot2::aes(x=pc_1, y=pc_2, color=sample)) +
      ggplot2::geom_point(size=point_size, alpha=point_alpha) +
      ggplot2::labs(x=base::paste0('PC 1: ', pca_variance[1], '%'),
                    y=base::paste0('PC 2: ', pca_variance[2], '%'),
                    title=title) +
      create_ggplot_theme() +
      ggplot2::scale_x_continuous(breaks=NULL) +
      ggplot2::scale_y_continuous(breaks=NULL) +
      ggplot2::guides(color=ggplot2::guide_legend(title=NULL,
                                                  label.theme=ggplot2::element_text(size=4),
                                                  keywidth=0.5,
                                                  keyheight=0.5))

  ggplot2::ggsave(plot=plot, out_path, width=width, height=height)
}

main <- function(in_file_name, out_file_name) {
  data_frame <- utils::read.table(file=in_file_name, sep='\t', header=TRUE)
  num_columns <- base::ncol(data_frame)
  ## first 2 columns are event_type, id
  sample_names <- base::colnames(data_frame)[3:num_columns]
  num_samples <- base::length(sample_names)
  psi_values <- data_frame[, 3:num_columns]
  pca_data <- calculate_pca(psi_values)
  pc_1 <- pca_data$pca$x[, 1]
  pc_2 <- pca_data$pca$x[, 2]
  plot_df <- base::data.frame(pc_1=pc_1, pc_2=pc_2, sample=sample_names)
  plot_pc_1_2(plot_df, pca_data$percent_variance, out_file_name)
}

main(in_file_name, out_file_name)
