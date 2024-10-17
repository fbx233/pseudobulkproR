#' Perfect Volcano Plot
#'
#' This function creates a customizable volcano plot for visualizing differentially expressed genes.
#'
#' @param plot_data A data frame containing the differential expression analysis results.
#' @param up_genes A data frame containing upregulated genes to be highlighted.
#' @param down_genes A data frame containing downregulated genes to be highlighted.
#' @param sig_genes A data frame containing significant genes to be labeled.
#' @param base_point_size Numeric value for the size of regular points. Default is 1.
#' @param sig_point_size Numeric value for the size of significant points. Default is 2.
#' @param point_alpha Numeric value for point transparency. Default is 0.6.
#' @param up_color Character string specifying color for upregulated genes. Default is "#ffad73".
#' @param down_color Character string specifying color for downregulated genes. Default is "#26b3ff".
#' @param ns_color Character string specifying color for non-significant genes. Default is "grey".
#' @param sig_up_fill Character string specifying fill color for significant upregulated points. Default is "red".
#' @param sig_down_fill Character string specifying fill color for significant downregulated points. Default is "steelblue".
#' @param line_type Character string specifying the type of threshold lines. Default is "dashed".
#' @param label_force Numeric value for label repulsion force. Default is 2.
#' @param label_nudge_y Numeric value for label Y-axis adjustment. Default is 1.
#' @param x_limits Numeric vector specifying x-axis limits. Default is c(-10, 10).
#' @param x_breaks Numeric vector specifying x-axis breaks. Default is seq(-10, 10, 2).
#' @param legend_pos Vector specifying legend position. Default is c(0.88, 0.89).
#' @param title_size Numeric value for title font size. Default is 10.
#' @param text_size Numeric value for text font size. Default is 9.
#' @param custom_labels Character vector of three elements for legend labels. Default uses gene counts.
#'
#' @return A ggplot object containing the volcano plot
#'
#' @example
#' perfect_volcano_plot(
#'   plot_data = deg_results,
#'   up_genes = up_deg,
#'   down_genes = down_deg,
#'   sig_genes = sig_deg,
#'   base_point_size = 1.2,
#'   up_color = "orange",
#'   down_color = "blue",
#'   x_limits = c(-8, 8)
#' )
#'@dataform
#' head(plot_data)
#'  symbol  fold_change adj_p_val gene_type
#'  1  Hspb3  0.001693839 0.0762149        ns
#'  2   Pbsn  1.048034935 1.0000000        ns
#'  3  Cdc45  0.772010765 0.5882763        ns
#' 
perfect_volcano_plot <- function(
  plot_data,
  up_genes,
  down_genes,
  sig_genes,
  base_point_size = 1,
  sig_point_size = 2,
  point_alpha = 0.6,
  up_color = "#ffad73",
  down_color = "#26b3ff",
  ns_color = "grey",
  sig_up_fill = "red",
  sig_down_fill = "steelblue",
  line_type = "dashed",
  label_force = 2,
  label_nudge_y = 1,
  x_limits = c(-10, 10),
  x_breaks = seq(-10, 10, 2),
  legend_pos = c(0.88, 0.89),
  title_size = 10,
  text_size = 9,
  custom_labels = NULL
) {
  
  # Generate default labels if not provided
  if(is.null(custom_labels)) {
      gene_counts <- table(plot_data$gene_type)
      custom_labels <- paste0(names(gene_counts), " ", gene_counts)
  }
  
  # Create the plot
  plot_data %>% 
      ggplot(aes(x = log2(fold_change), y = -log10(adj_p_val))) + 
      geom_point(aes(color = gene_type), 
                alpha = point_alpha, 
                shape = 16, 
                size = base_point_size) + 
      geom_point(data = up_genes, 
                shape = 21, 
                size = sig_point_size, 
                fill = sig_up_fill, 
                colour = "black") + 
      geom_point(data = down_genes, 
                shape = 21, 
                size = sig_point_size, 
                fill = sig_down_fill, 
                colour = "black") + 
      geom_hline(yintercept = -log10(0.05), 
                linetype = line_type) + 
      geom_vline(xintercept = c(log2(0.5), log2(2)), 
                linetype = line_type) +
      geom_label_repel(data = sig_genes, 
                      aes(label = symbol), 
                      force = label_force, 
                      nudge_y = label_nudge_y, 
                      max.overlaps = Inf) +
      scale_color_manual(values = c("up" = up_color, 
                                  "down" = down_color, 
                                  "ns" = ns_color),
                       labels = custom_labels) + 
      scale_x_continuous(breaks = x_breaks, 
                       limits = x_limits) +
      labs(x = "log2(fold change)",
           y = "-log10(adjusted P-value)", 
           colour = "Expression change") +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      theme_bw() +   
      theme(
          panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 0.5),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.title = element_text(face = "bold", 
                                  color = "black", 
                                  size = title_size),
          axis.text = element_text(color = "black", 
                                 size = text_size, 
                                 face = "bold"),
          legend.background = element_blank(),
          legend.title = element_text(face = "bold", 
                                    color = "black", 
                                    size = title_size),
          legend.text = element_text(face = "bold", 
                                   color = "black", 
                                   size = text_size),
          legend.spacing.x = unit(0, "cm"),
          legend.position = "inside",
          legend.position.inside = legend_pos,
          legend.justification = c(0.5, 0.5)
      )
}