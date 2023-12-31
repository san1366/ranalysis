readNanoStringRccSet <- function(rccFiles, rlfFile = NULL, phenoDataFile = NULL, phenoDataRccColName = "^RCC", 
    phenoDataColPrefix = "") {
    data <- structure(lapply(rccFiles, readRccFile), names = basename(rccFiles))
    assay <- do.call(cbind, lapply(data, function(x) structure(x[["Code_Summary"]][["Count"]], 
        names = rownames(x[["Code_Summary"]]))))
    if (is.null(phenoDataFile)) {
        pheno <- annotatedDataFrameFrom(assay, byrow = FALSE)
    }
    else {
        pheno <- read.csv(phenoDataFile, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
        j <- grep(phenoDataRccColName, colnames(pheno), ignore.case = TRUE)
        if (length(j) == 0L) 
            stop("Column `phenoDataRccColName` not found in `phenoDataFile`")
        else if (length(j) > 1L) 
            stop("Multiple columns in `phenoDataFile` match `phenoDataRccColName`")
        missingPhenoCount <- sum(!(colnames(assay) %in% pheno[[j]]))
        rownames(pheno) <- pheno[[j]]
        pheno[[j]] <- NULL
        pheno <- pheno[colnames(assay), , drop = FALSE]
        if (missingPhenoCount != 0L) {
            rownames(pheno) <- colnames(assay)
            warning(sprintf("Column `phenoDataRccColName` in `phenoDataFile` is missing %d of %d Samples", 
                missingPhenoCount, ncol(assay)))
        }
        if (phenoDataColPrefix != "") {
            colnames(pheno) <- paste0(phenoDataColPrefix, colnames(pheno))
        }
        pheno <- AnnotatedDataFrame(pheno, dimLabels = c("sampleNames", "sampleColumns"))
    }
    feature <- lapply(data, function(x) {
        x[["Code_Summary"]][, c("CodeClass", "GeneName", "Accession")]
    })
    stopifnot(all(vapply(feature, function(x) identical(feature[[1L]], x), FUN.VALUE=logical(1))))
    feature <- feature[[1L]]
    feature[["IsControl"]] <- NA
    feature[["IsControl"]][feature[["CodeClass"]] %in% .codeClassMetadata[.codeClassMetadata[["IsControl"]], 
        "CodeClass"]] <- TRUE
    feature[["IsControl"]][feature[["CodeClass"]] %in% .codeClassMetadata[!.codeClassMetadata[["IsControl"]], 
        "CodeClass"]] <- FALSE
    feature[["ControlConc"]] <- NA_real_
    feature[["ControlConc"]][feature[["CodeClass"]] == "Negative"] <- 0
    feature[["ControlConc"]][feature[["CodeClass"]] == "Positive"] <- as.numeric(sub("POS_\\w\\((.*)\\)", 
        "\\1", feature[["GeneName"]][feature[["CodeClass"]] == "Positive"]))
    if (is.null(rlfFile)) {
        rlfHeader <- list()
    }
    else if (!is.null(rlfFile)) {
        rlfData <- readRlfFile(rlfFile)
        rlfHeader <- metadata(rlfData)
        rlfHeader[["RlfFileDate"]] <- as.character(rlfHeader[["RlfFileDate"]])
        rlfData <- as.data.frame(rlfData)
        rlfData <- rlfData[rlfData[["CodeClassActive"]] %in% c(2L, 3L), , drop = FALSE]
        rownames(rlfData) <- sprintf("%s_%s_%s", rlfData[["CodeClass"]], rlfData[["GeneName"]], 
            rlfData[["Accession"]])
        if (!identical(sort(rownames(feature)), sort(rownames(rlfData)))) 
            stop("featureData mismatch between RLF and RCC files")
        rlfData <- rlfData[rownames(feature), , drop = FALSE]
        for (j in c("CodeClass", "GeneName", "Accession", "CodeClassActive")) {
            rlfData[[j]] <- NULL
        }
        feature <- cbind(feature, rlfData)
    }
    feature <- AnnotatedDataFrame(feature, dimLabels = c("featureNames", "featureColumns"))
    name <- unlist(lapply(data, function(x) x[["Sample_Attributes"]][["SampleOwner"]]))
    name <- unique(na.omit(name))
    experiment <- MIAME(name = name, other = rlfHeader)
    annotation <- unlist(lapply(data, function(x) x[["Sample_Attributes"]][["GeneRLF"]]))
    annotation <- unique(annotation)
    if (length(annotation) > 1L) 
        stop("RCC files do not have the same GeneRLF attribute")
    protocol <- do.call(rbind, lapply(seq_along(rccFiles), function(i) {
        x <- data[[i]][["Sample_Attributes"]]
        x <- x[, setdiff(names(x), "GeneRLF")]
        cbind(data[[i]][["Header"]], x, data[[i]][["Lane_Attributes"]])
    }))
    protocol <- AnnotatedDataFrame(protocol, .rccMetadata[["protocolData"]], dimLabels = c("sampleNames", 
        "sampleColumns"))
    NanoStringRccSet(assayData = assay, phenoData = pheno, featureData = feature, experimentData = experiment, 
        annotation = annotation, protocolData = protocol, check = FALSE)
}