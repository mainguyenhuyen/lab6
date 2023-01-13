#!/usr/bin/env nextflow

params.logFC = 1
params.FDR = 0.05

process differentialAnalysis {
  input:
    path countFile
  
  output:
    path 'analysis_results.xlsx'

  script:
    """
      #!/usr/bin/env R -e "\
      library(edgeR); \
      library(openxlsx); \
      countTable <- read.table('$countFile', header = TRUE, as.is = TRUE, row.names = 1, sep = '\\t'); \
      dge <- DGEList(countTable); \
      dge <- calcNormFactors(dge); \
      treatment = factor(rep(c('CHIR', 'CTRL'), each = 3)); \
      pairs <- factor(rep(c('A', 'B', 'C'), 2)); \
      designDF <- data.frame(Treatment = treatment, Pair = pairs); \
      design <- model.matrix(~ 0 + Treatment + Pair, data = designDF);  \
      dge <- estimateDisp(dge, design, robust = TRUE); \
      fit <- glmQLFit(dge, design, robust = TRUE); \
      contrastMat <- makeContrasts(TreatmentVsControl = TreatmentCHIR - TreatmentCTRL, levels = design); \
      qlfRes <- glmQLFTest(fit, contrast = contrastMat); \
      topRes <- topTags(qlfRes, n = nrow(fit[['counts']])); \
      topRes <- subset(topRes[['table']], abs(logFC) > ${params.logFC} & FDR < ${params.FDR}); \
      write.xlsx(topRes, 'analysis_results.xlsx', rowNames = TRUE); "
    """
}

workflow {
  inputFile = Channel.fromPath(params.input)
  differentialAnalysis(inputFile)
}
