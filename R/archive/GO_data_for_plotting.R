#'generate region metadata from regions of interest given provided annotation files
#'
#' @param candidate_list genes_of_interest
#' @param bg_genes background genes
#' @param db database
#' @param go_ids gene ids to be used
#' @return Function returns dataframe of enriched gene ontology categories with their p-values
#' @export GO_data_for_plotting

GO_data_for_plotting <- function(candidate_list,
                                 bg_genes,
                                 db= useMart('ENSEMBL_MART_ENSEMBL',dataset='clfamiliaris_gene_ensembl', host="https://www.ensembl.org"),
                                 go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'),
                                               filters='external_gene_name', values=bg_genes, mart=db)
) {

  bg_genes=pqlseq_res_oi$gene
  # Read in genes of interest'
  candidate_list <- pqlseq_res_oi_fdr1$gene[pqlseq_res_oi_fdr1$beta > 0]

  # build the gene 2 GO annotation list (needed to create topGO object)
  gene_2_GO=unstack(go_ids[,c(1,2)])

  # remove any candidate genes without GO annotation
  keep = candidate_list %in% go_ids[,2]
  keep =which(keep==TRUE)
  candidate_list=candidate_list[keep]

  # make named factor showing which genes are of interest
  geneList=factor(as.integer(bg_genes %in% candidate_list))
  names(geneList)= bg_genes

  GOdata=new('topGOdata', ontology='BP', allGenes = geneList,
             annot = annFUN.gene2GO, gene2GO = gene_2_GO)

  classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher')

  # define test using the weight01 algorithm (default) with fisher
  weight_fisher_result=runTest(GOdata, algorithm='weight01', statistic='fisher')

  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO=usedGO(GOdata)
  all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

  #performing BH correction on our p values
  p.adj=p.adjust(all_res$weightFisher,method="BH")

  # create the file with all the statistics from GO analysis
  all_res_final=cbind(all_res,p.adj)
  all_res_final=all_res_final[order(all_res_final$p.adj),]

  #graphing the results
  all_res_final_ggplot <- all_res_final
  all_res_final_ggplot$GO.id_labeled <- paste0(all_res_final_ggplot[,2], " (", all_res_final_ggplot[,1], ")")

  all_res_final_ggplot$GO.id_labeled <- factor(all_res_final_ggplot$GO.id_labeled,                                    # Factor levels in decreasing order
                                               levels = all_res_final_ggplot$GO.id_labeled[order(all_res_final_ggplot$p.adj, decreasing = TRUE)])

  return(all_res_final_ggplot)
}
