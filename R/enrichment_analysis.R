#' @importFrom AnnotationDbi select
#' @importFrom clusterProfiler enrichGO enrichKEGG enrichWP
#' @importFrom DOSE setReadable
#' @importFrom org.Dr.eg.db org.Dr.eg.db
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom stringr str_split
#'
NULL

#' Perform enrichment analysis on a set of Entrez IDs
#'
#' This function performs enrichment analysis on a set of Entrez IDs.
#'
#' @param entrezIDs A character vector of gene Entrez IDs.
#' @inheritParams entrezGenes
#' @param funString Name of enrichment function from \code{clusterProfiler}.
#' Must be a character selected from 'enrichGO', 'enrichKEGG' and 'enrichWP'.
#'
#' @return Enrichment result.
#'
#' @keywords internal
#'
getEnrichmentResult <- function(entrezIDs,
                                species = c('human', 'mouse', 'zebrafish'),
                                funString = c('enrichGO',
                                              'enrichKEGG',
                                              'enrichWP')){
    species <- match.arg(species, c('human', 'mouse', 'zebrafish'))
    funString <- match.arg(funString, c('enrichGO', 'enrichKEGG', 'enrichWP'))
    df <- data.frame(database = c('org.Hs.eg.db', 'org.Mm.eg.db',
                                  'org.Dr.eg.db'),
                     code = c('hsa', 'dre', 'mmu'))
    rownames(df) <- c('human', 'mouse', 'zebrafish')
    db <- df[species, 1]
    code <- df[species, 2]
    return(switch(funString,
           "enrichGO" = enrichGO(entrezIDs, OrgDb=db, ont="BP",
                                 readable=TRUE, pvalueCutoff=0),
           "enrichWP" = setReadable(enrichWP(entrezIDs, 'Homo sapiens'),
                                    db, 'ENTREZID'),
           "enrichKEGG" = setReadable(enrichKEGG(entrezIDs, code),
                                      db, 'ENTREZID'),
    ))
}

#' Convert gene symbols to Entrez IDs.
#'
#' This function converts gene symbols to Entrez IDs.
#'
#' @param genes A character vector of gene symbols.
#' @param species Species. Must be one of 'human', 'mouse' and 'zebrafish'.
#'
#' @return The Entrez IDs of the genes.
#'
#' @keywords internal
#'
entrezGenes <- function(genes, species = c('human', 'mouse', 'zebrafish')){
    species <- match.arg(species, c('human', 'mouse', 'zebrafish'))
    return(switch(species,
           'human' = AnnotationDbi::select(org.Hs.eg.db,
                                           keys=genes,
                                           columns=c("ENTREZID", "SYMBOL"),
                                           keytype="SYMBOL")[[2]],
           'zebrafish' = AnnotationDbi::select(org.Dr.eg.db,
                                               keys=genes,
                                               columns=c("ENTREZID", "SYMBOL"),
                                               keytype="SYMBOL")[[2]],
           'mouse' = AnnotationDbi::select(org.Mm.eg.db, keys=genes,
                                           columns = c("ENTREZID", "SYMBOL"),
                                           keytype="SYMBOL")[[2]]))
}


#' Perform enrichment analysis on a set of genes
#'
#' This function performs enrichment analysis on a set of genes.
#'
#' @inheritParams entrezGenes
#' @inheritParams getEnrichmentResult
#'
#' @return Enrichment result.
#'
#' @export
#'
genesER <- function(genes, species,
                    funString = c('enrichGO','enrichKEGG', 'enrichWP'))
    return(getEnrichmentResult(entrezGenes(genes, species), species, funString))

#' Extract genes enriched for term
#'
#' This function extracts genes enriched for term from an \code{enrichResult}
#' object.
#'
#' @param er Enrichment result.
#' @param terms Terms for which enriched genes should be extracted.
#' @param negTerms Terms for which enriched genes should be subtracted from the
#' genes enriched for \code{terms}.
#'
#' @return Genes enriched for terms.
#'
#' @export
#'
termGenes <- function(er, terms, negTerms = NULL){
    posGenes <- str_split(er@result[er@result$Description %in%
                                        terms, ]$geneID, '/')[[1]]
    if(!is.null(negTerms))
        negGenes <- str_split(er@result[er@result$Description %in%
                                            negTerms, ]$geneID, '/')[[1]] else
                                                negGenes <- NULL
    return(sort(setdiff(posGenes, negGenes)))
}
