#' @title barplot for enrichment result
#'
#' @description
#' plot barplot for enrichment result by use \code{ggplot2} package.
#'
#' @author Abao Xing \url{https://github.com/GitBioinformatics}
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param obj a enriched result in data.frame format.
#' @param x variable for x-axis, one of 'GeneRatio' (default), 'BgRatio' and 'count.enriched'.
#' @param color.by variable that used to color enriched terms,
#' 'arg' should be one of 'pvalue', 'p.adjust' (default), 'qvalue'.
#' @param order.by The order of the Y-axis, NULL (default) means no sorting.
#' @param order.decreasing logical. Should the order.by order be increasing or decreasing?
#' This param is invalid when order.by = NULL.
#' @param order.show A number or a vector of terms. If it is a number,
#' the first n terms will be displayed. If it is a vector of terms,
#' the selected terms will be displayed.
#' @param fill.color.range a vector of two color values, set the colors.
#' @param font.sizefont size
#' @param x.title A character string or expression indicating a title of x-axis If NULL, the title is shown by default.
#' @param y.title A character string or expression indicating a title of y-axis If NULL, the title is shown by default.
#' @param legend.size.title A character string or expression indicating a title of guide (size). If NULL, the title is shown by default.
#' @param legend.color.title A character string or expression indicating a title of guide (color). If NULL, the title is shown by default.
#' @param legend.position the position of legends ("none", "left", "right" (default), "bottom", "top", or two-element numeric vector [0, 1])
#' @param legend.direction layout of items in legends, one of "horizontal" or "vertical" (default).
#'
#' @return a \code{ggplot2} instance
#'
#' @examples
#' # ---------------------------------------------------------------------------
#' # example one
#' library(OmicsPEA)
#'
#' data(Metabolites)
#' data(Metabolites.Pathways)
#' colnames(Metabolites.Pathways) <- c('Metabolites', 'itemId', 'pathwayName')
#' PEA.RES = PEA(omics.id.set = Metabolites$ID, omics.pathway.db = Metabolites.Pathways, min.enrich.size = 1)
#'
#' barplot(
#'   obj = PEA.RES,
#'   x = 'omicsRatio',
#'   color.by = 'p.adjust',
#'   order.by = 'x',
#'   order.decreasing = FALSE,
#'   order.show = c('Phospholipid Metabolism',
#'                  'Phospholipid Metabolism',
#'                  'Glycosphingolipid Metabolism',
#'                  'Vitamin B3 (Niacin, NAD+) Metabolism',
#'                  'Krebs Cycle',
#'                  'Histidine, Histamine, Carnosine Metabolism')
#' )
#'
#' # ---------------------------------------------------------------------------
#' # example two
#' library(OmicsPEA)
#' library(org.Hs.eg.db)
#' library(AnnotationDbi)
#'
#' data(Gene.Expression.Demo)
#' data(KEGG.Pathways)
#' colnames(KEGG.Pathways) <- c('itemId', 'gene_symbol', 'PathwayID', 'pathwayName')
#' Gene.Expression.Demo$gene_id <- mapIds(org.Hs.eg.db, Gene.Expression.Demo$gene_symbol, column = "ENTREZID", keytype = "SYMBOL")
#' PEA.RES = PEA(omics.id.set = Gene.Expression.Demo$gene_id, omics.pathway.db = KEGG.Pathways, min.enrich.size = 10, ratio.numeric = FALSE)
#'
#' colnames(PEA.RES)
#' barplot(
#'   obj = PEA.RES,
#'   x = "count.enriched",
#'   color.by = 'p.adjust',
#'   order.by = 'x',
#'   order.decreasing = TRUE,
#'   order.show = 10,
#'   font.size = 10,
#'   x.title = 'Enriched Count',
#'   y.title = NULL,
#'   legend.size.title = NULL,
#'   legend.color.title = NULL,
#'   legend.position = 'right',
#'   legend.direction = 'vertical'
#' )
#'
barplot <- function(
    obj,
    x = "omicsRatio",
    color.by = 'p.adjust',
    order.by = NULL,
    order.decreasing = FALSE,
    order.show = 10,
    fill.color.range = c('red', 'blue'),
    font.size = 12,
    x.title = NULL,
    y.title = NULL,
    legend.size.title = NULL,
    legend.color.title = NULL,
    legend.position = 'right',
    legend.direction = 'vertical'
) {

  obj[['omicsRatio']] <- sapply(obj[['omicsRatio']], function(x) eval(parse(text = x)))
  obj[['BgRatio']] <- sapply(obj[['BgRatio']], function(x) eval(parse(text = x)))

  color.by <- match.arg(color.by, c('pvalue', 'p.adjust', 'qvalue'))
  x <- match.arg(x, c('omicsRatio', 'BgRatio', 'count.enriched'))

  legend.direction <- match.arg(legend.direction, c('horizontal', 'vertical'))
  if (class(legend.position) == 'character' && length(legend.position) == 1) {
    legend.position <- match.arg(legend.position, c('none', 'left', 'bottom', 'right', 'top'))
  } else if (class(legend.position) == 'numeric') {
    if (length(legend.position) != 2) {
      stop('the input of legend.position is wrong, the value can be one of "none", "left", "right" (default), "bottom", "top", or two-element numeric vector [0, 1]')
    } else {
      if (legend.position[1] < 0 || legend.position[1] > 1 || legend.position[2] < 0 || legend.position[2] > 1) {
        stop('the input of legend.position is wrong, the value can be one of "none", "left", "right" (default), "bottom", "top", or two-element numeric vector [0, 1]')
      }
    }
  } else {
    stop('the input of legend.position is wrong, the value can be one of "none", "left", "right" (default), "bottom", "top", or two-element numeric vector [0, 1]')
  }

  if (!is.null(order.by)) {
    if (order.by !=  'x' && !order.by %in% colnames(obj)) { message('wrong "order.by" parameter; set to default `order.by = "x"`'); order.by <- "x" }
    if (order.by == "x") { obj <- dplyr::mutate(obj, x = eval(parse(text = x))) }
    idx <- order(obj[[order.by]], decreasing = order.decreasing)
    obj$ID <- factor(obj$ID, levels = rev(unique(obj$ID[idx])))
  }

  if (class(order.show) == 'numeric' && length(order.show) == 1) {
    obj <- obj %>% slice(1:order.show)
  } else if(length(order.show) > 1) {
    obj <- obj %>% dplyr::filter(ID %in% order.show)
  }

  gg <- ggplot(data = obj) +
    geom_bar(aes_string(x = 'ID', y = x, fill = color.by), stat = 'identity') +
    coord_flip() +
    scale_fill_gradient(low = fill.color.range[1], high = fill.color.range[2]) +
    xlab(ifelse(test = is.null(y.title), yes = 'Pathway Name', no = y.title)) +
    ylab(ifelse(test = is.null(x.title), yes = x, no = x.title)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0))

  gg <- gg + theme_bw() + theme(
    legend.title = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 0, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(t = 0, r = 0, b = 5.5, l = 0, unit = 'pt')),
    legend.text = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 0, vjust = 0.5, angle = 0.5, lineheight = 0.9, margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'pt')),
    axis.text.y.left = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 1, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(t = 0, r = 5.5, b = 0, l = 0, unit = 'pt')),
    axis.text.x.bottom = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 0.5, vjust = 0, angle = 0, lineheight = 0.9, margin = margin(t = 5.5, r = 0, b = 0, l = 0, unit = 'pt')),
    axis.title.y.left = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 0.5, vjust = 0.5, angle = 90, lineheight = 0.9, margin = margin(t = 5.5, r = 0, b = 0, l = 0, unit = 'pt')),
    axis.title.x.bottom = element_text(family = NA, face = 'plain', colour = '#000000', size = font.size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9, margin = margin(t = 5.5, r = 0, b = 0, l = 0, unit = 'pt')),
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.justification = "right"
  )

  gg
}
