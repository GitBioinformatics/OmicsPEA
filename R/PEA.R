#' @title Pathway Enrichment Analysis (PEA)
#'
#' @description
#' PEA (Pathway Enrichment Analysis) of a gene/protein/metabolite or other set.
#' This function will do PEA and return the enriched results.
#'
#' @param omics.id.set required, a vector of gene/protein/metabolite or other omics id.
#' @param omics.pathway.db required, It is a pathway database in data.frame format and must contain two columns, \code{itemId} and \code{pathwayName}.
#' @param adjust.p.method character, one of "none", "BH (default)", "bonferroni", "fdr", "hochberg", "holm", "hommel", "BY", select method used in p.adjust.
#' @param min.enrich.size numeric, minimal size of \code{omics.id.set} annotated by \code{omics.pathway.db}, 10 is default.
#' @param max.enrich.size numeric, maximal size of \code{omics.id.set} annotated by \code{omics.pathway.db}, 1000 is default.
#' @param ratio.numeric logical, ratio value is in character or numeric format, defaults to FALSE (character).
#'
#' @return a \code{enrich.res} data.frame instance
#'
#' @author Abao Xing \url{https://github.com/GitBioinformatics}
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import tibble
#'
#' @examples
#' # ---------------------------------------------------------------------------
#' # example one
#' library(OmicsPEA)
#'
#' data(Metabolites)
#' data(Metabolites.Pathways)
#' colnames(Metabolites.Pathways) <- c('Metabolites', 'itemId', 'pathwayName')
#' enrich.res = PEA(omics.id.set = Metabolites$ID, omics.pathway.db = Metabolites.Pathways, min.enrich.size = 1)
#' head(enrich.res)
#'
#' # ---------------------------------------------------------------------------
#' # example two
#' # use gene_symbol as itemId
#' library(OmicsPEA)
#'
#' data(Gene.Expression.Demo)
#' data(Wiki.Pathways)
#' colnames(Wiki.Pathways) <- c('GeneID', 'itemId', 'PathwayID', 'pathwayName')
#' enrich.res = PEA(omics.id.set = Gene.Expression.Demo$gene_symbol, omics.pathway.db = Wiki.Pathways, min.enrich.size = 10)
#' head(enrich.res)
#'
#' # ---------------------------------------------------------------------------
#' # example three
#' # use gene_id as itemId
#' library(OmicsPEA)
#' library(org.Hs.eg.db)
#' library(AnnotationDbi)
#'
#' data(Gene.Expression.Demo)
#' data(KEGG.Pathways)
#' colnames(KEGG.Pathways) <- c('itemId', 'gene_symbol', 'PathwayID', 'pathwayName')
#' Gene.Expression.Demo$gene_id <- mapIds(org.Hs.eg.db, Gene.Expression.Demo$gene_symbol, column = "ENTREZID", keytype = "SYMBOL")
#' enrich.res = PEA(omics.id.set = Gene.Expression.Demo$gene_id, omics.pathway.db = KEGG.Pathways, min.enrich.size = 10)
#' head(enrich.res)
#'
#' @export
PEA <- function (
    omics.id.set,
    omics.pathway.db,
    adjust.p.method = "BH",
    min.enrich.size = 10,
    max.enrich.size = 1000,
    ratio.numeric = FALSE
) {

  # check omics.id.set is or not in omics.pathway.db
  omics.id.set <- unique(as.character(omics.id.set))
  omics.id.set <- intersect(omics.id.set, as.character(omics.pathway.db$itemId))

  # Convert omics.pathway.db to Omics2Pathway for format
  Omics2Pathway <- omics.pathway.db %>%
    group_by(itemId) %>%
    summarize(pathwayName = list(pathwayName)) %>%
    ungroup() %>%
    mutate(pathwayName = map(pathwayName, ~unique(.x))) %>%
    deframe()

  # Convert omics.pathway.db to Pathway2Omics for format
  Pathway2Omics <- omics.pathway.db %>%
    group_by(pathwayName) %>%
    summarize(itemId = list(itemId)) %>%
    ungroup() %>%
    mutate(itemId = map(itemId, ~unique(.x))) %>%
    deframe()

  omics.id.set <- as.character(unique(omics.id.set))
  qPathwayID2OmicsID <- Omics2Pathway[omics.id.set]
  qOmicsID <- unlist(qPathwayID2OmicsID)
  qPathwayID2OmicsID.df <- data.frame(pathwayID = rep(names(qPathwayID2OmicsID), times = lapply(qPathwayID2OmicsID, length)), omicsID = qOmicsID)
  qPathwayID2OmicsID.df <- unique(qPathwayID2OmicsID.df)
  qOmicsID2PathwayID <- with(qPathwayID2OmicsID.df, split(as.character(pathwayID), as.character(omicsID)))
  pathwayID <- names(Omics2Pathway)

  qOmicsID2PathwayID <- lapply(qOmicsID2PathwayID, intersect, pathwayID)
  qOmicsID <- unique(names(qOmicsID2PathwayID))
  tmp.list <- lapply(Pathway2Omics, function(x) length(intersect(x, omics.id.set)) > 0)
  omicsID2PathwayID <- Pathway2Omics[unlist(tmp.list)]
  omicsID2PathwayID <- lapply(omicsID2PathwayID, intersect, pathwayID)
  geneSets <- omicsID2PathwayID
  idx <- DOSE:::get_geneSet_index(omicsID2PathwayID, min.enrich.size, max.enrich.size)

  omicsID2PathwayID <- omicsID2PathwayID[idx]
  qOmicsID2PathwayID <- qOmicsID2PathwayID[idx]
  qOmicsID <- unique(names(qOmicsID2PathwayID))
  k <- sapply(qOmicsID2PathwayID, length)
  k <- k[qOmicsID]
  M <- sapply(omicsID2PathwayID, length)
  M <- M[qOmicsID]
  N <- rep(length(pathwayID), length(M))
  n <- rep(length(qPathwayID2OmicsID), length(M))
  args.df <- data.frame(numWdrawn = k - 1, numW = M, numB = N - M, numDrawn = n)
  pvalues <- apply(args.df, 1, function(n) phyper(n[1], n[2], n[3], n[4], lower.tail = FALSE))
  omicsRatio <- apply(data.frame(a = k, b = n), 1, function(x) paste(x[1], "/", x[2], sep = "", collapse = ""))
  BgRatio <- apply(data.frame(a = M, b = N), 1, function(x) paste(x[1], "/", x[2], sep = "", collapse = ""))
  Over <- data.frame(ID = as.character(qOmicsID), omicsRatio = omicsRatio, BgRatio = BgRatio, pvalue = pvalues, stringsAsFactors = FALSE)
  p.adj <- p.adjust(Over$pvalue, method = adjust.p.method)
  qobj <- tryCatch(qvalue::qvalue(p = Over$pvalue, lambda = 0.05, pi0.method = "bootstrap"), error = function(e) NULL)
  if (inherits(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  omicsID <- sapply(qOmicsID2PathwayID, function(i) paste(i, collapse = "/"))
  omicsID <- omicsID[qOmicsID]
  Over <- data.frame(Over, p.adjust = p.adj, qvalue = qvalues, omicsID = omicsID, count.enriched = k, stringsAsFactors = FALSE)

  nc <- ncol(Over)
  Over <- Over[, c(1, nc, 2:(nc - 1))]
  Over <- Over[order(pvalues), ]
  Over$ID <- as.character(Over$ID)
  row.names(Over) <- as.character(Over$ID)

  if (ratio.numeric) {
    Over[['omicsRatio']] <- sapply(Over[['omicsRatio']], function(x) eval(parse(text = x)))
    Over[['BgRatio']] <- sapply(Over[['BgRatio']], function(x) eval(parse(text = x)))
  }
  return(Over)
}
