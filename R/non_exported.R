#'@importFrom stats na.omit lm fitted predict median
#'@importFrom matrixStats colMedians
#'@import drc
standard_names <- function(chans, reps, dataType) {
    if (dataType == "intensity") {
        
        tempNames <- paste("rep", rep(seq_len(reps), each = chans), "_", rep(paste("C", c("ontrol", 
            0:(chans - 2)), sep = ""), reps), sep = "")
        
        
        
    } else {
        tempNames <- paste("rep", rep(seq_len(reps), each = chans), "_", rep(paste("C", 0:(chans - 
            1), sep = ""), reps), sep = "")
        
    }
    return(tempNames)
}


index_matrix <- function(nchans, reps, dataType, dataFrame) {
    
    
    if (reps <= 1) {
        NULL
    } else {
        
        if (dataType == "intensity") {
            channe <- seq_len(((nchans - 1) * reps))
            chans <- nchans - 1
            ser <- standard_names(nchans, reps, dataType)
            ser <- ser[-seq(1, by = nchans, length.out = reps)]
            
        } else {
            ser <- standard_names(nchans, reps, dataType)
            channe <- seq_len((nchans * reps))
            chans <- nchans
            
        }
        
        combinations <- t(utils::combn(reps, 2))
        combmat <- matrix(rep(as.vector(t(combinations)), chans), ncol = 2, byrow = TRUE)
        
        # create factor for repeat ...
        
        name.vec <- seq_len((chans * reps))
        repfac <- rep(seq_len(chans), times = reps)
        index <- rep(0:(reps - 1), each = chans)
        combfac <- rep(seq_len(reps), each = chans)
        
        columnindex <- matrix(0, ncol = 5, nrow = chans * nrow(combinations))
        
        colnames(columnindex) <- c("concentration", "rep1", "rep2", "index1", "index2")
        columnindex[, 1] <- rep(seq_len(chans), each = nrow(combinations))
        
        columnindex[, 2:3] <- combmat
        
        columnindex
        # create matrix which will be indexed by first 3 columns of column index
        
        index.mat <- matrix(name.vec, ncol = reps)
        
        columnindex[, 4] <- index.mat[columnindex[, c(1, 2)]]
        columnindex[, 5] <- index.mat[columnindex[, c(1, 3)]]
        
        create.names <- rep("", nrow(columnindex))
        
        
        create.names <- paste(ser[columnindex[, 4]], "vs", ser[columnindex[, 5]])
        
        final.mat <- data.frame(names = create.names, columnindex)
        final.mat$concentration <- final.mat$concentration - 1
        final.mat[!is.na(rowSums(matrix(match(columnindex, channe), ncol = 5))), ]
    }
    
}


peptide_match <- function(dr1, dr2, nchan) {
    maxrow <- max(nrow(dr1), nrow(dr2))
    minrow <- min(nrow(dr1), nrow(dr2))
    
    adVal <- maxrow - minrow
    
    if (nrow(dr1) == maxrow) {
        dr2$addedVals <- adVal
        big.pep <- dr1
        small.pep <- dr2
    } else {
        dr1$addedVals <- adVal
        big.pep <- dr2
        small.pep <- dr1
    }
    
    
    newframe <- big.pep
    colnames(newframe) <- colnames(small.pep)
    newframe[seq_len(minrow), ] <- small.pep
    
    
    newframe[(minrow + 1):maxrow, ] <- base::colMeans(small.pep)
    
    if (all.equal(dim(big.pep), dim(dr1)) == TRUE) {
        
        dr2 <- newframe
    } else {
        dr1 <- newframe
    }
    
    list(dr1 = dr1, dr2 = dr2)
}


uniprotGene <- function(organism) {
    
    
    if (organism == "H.sapiens") {
        query = "<query model=\"genomic\" view=\"Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol\" sortOrder=\"Protein.primaryAccession ASC\" >
    <constraint path=\"Protein.organism.shortName\" op=\"=\" value=\"H. sapiens\" code=\"A\" />
    </query>"
        
        ret = httr::POST("http://www.humanmine.org/humanmine/service/query/results", body = list(query = query, 
            format = "json"), encode = "form")
    } else if (organism == "D.melanogaster") {
        
        query = "<query model=\"genomic\" view=\"Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol\" sortOrder=\"Protein.primaryAccession ASC\" >
    <constraint path=\"Protein.organism.shortName\" op=\"=\" value=\"D. melanogaster\" code=\"A\" />
    </query>"
        
        ret = httr::POST("http://www.flymine.org/flymine/service/query/results", body = list(query = query, 
            format = "json"), encode = "form")
        
    } else if (organism == "M.musculus") {
        query = "<query model=\"genomic\" view=\"Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol\" sortOrder=\"Protein.primaryAccession ASC\" >
    <constraint path=\"Protein.organism.shortName\" op=\"=\" value=\"M. musculus\" code=\"A\" />
    </query>"
        ret = httr::POST("http://www.mousemine.org/mousemine/service/query/results", body = list(query = query, 
            format = "json"), encode = "form")
    } else if (organism == "R.norvegicus") {
        query = "<query model=\"genomic\" view=\"Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol\" sortOrder=\"Protein.primaryAccession ASC\" >
    <constraint path=\"Protein.organism.shortName\" op=\"=\" value=\"R. norvegicus\" code=\"A\" />
    </query>"
        
        ret = httr::POST("http://www.ratmine.org/ratmine/service/query/results", body = list(query = query, 
            format = "json"), encode = "form")
    } else {
        query = "<query model=\"genomic\" view=\"Protein.primaryAccession Protein.uniprotAccession Protein.genes.symbol\" sortOrder=\"Protein.primaryAccession ASC\" >
    <constraint path=\"Protein.organism.shortName\" op=\"=\" value=\"C. elegans\" code=\"A\" />
    </query>"
        
        ret = httr::POST("http://www.humanmine.org/humanmine/service/query/results", body = list(query = query, 
            format = "json"), encode = "form")
    }
    
    response = jsonlite::fromJSON(httr::content(ret, as = "text"))
    
    data.prots <- response$results[, c(1, 3)]
    colnames(data.prots) <- c("Entry", "Gene.names")
    as.data.frame(data.prots[, seq_len(2)], stringsAsFactors = FALSE)
    
}


panel_shadeNtext <- function(x, y, corr = NULL, col.regions, ...) {
    if (is.null(corr)) 
        corr <- stats::cor(x, y, use = "
                       pair")
    ncol <- 14
    pal <- col.regions(ncol)
    col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, length = ncol + 1), include.lowest = TRUE))
    usr <- graphics::par("usr")
    graphics::rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], border = NA)
    graphics::box(col = "lightgray")
    on.exit(graphics::par(usr))
    graphics::par(usr = c(0, 1, 0, 1))
    r <- formatC(corr, digits = 2, format = "f")
    cex.cor <- 0.8/graphics::strwidth("-X.xx")
    graphics::text(0.5, 0.5, r, cex = cex.cor)
}

shape_for_ggplot_pred <- function(df_ordered, conc, pred.names) {
    cols_to_keep_pred <- c(pred.names, "GeneID", "Accession")
    
    forggplot_pred <- vector(mode = "list", length = length(df_ordered$GeneID))
    
    for (i in seq_len(length(df_ordered$GeneID))) {
        tmp_pred <- df_ordered[, cols_to_keep_pred]
        forggplot_pred[[i]] <- reshape2::melt(tmp_pred[i, ], id = c("GeneID", "Accession"), na.rm = TRUE)
    }
    
    forggplot_pred_1 <- do.call(rbind, forggplot_pred)
    forggplot_pred_1 <- data.frame(forggplot_pred_1, x = conc)
}

shape_for_ggplot_perc <- function(df_ordered, conc, finalNames) {
    cols_to_keep_perc <- c(finalNames, "GeneID", "Accession")
    
    forggplot_perc <- vector(mode = "list", length = length(df_ordered$GeneID))
    for (i in seq_len(length(df_ordered$GeneID))) {
        tmp_perc <- df_ordered[, cols_to_keep_perc]
        forggplot_perc[[i]] <- reshape2::melt(tmp_perc[i, ], id = c("GeneID", "Accession"), na.rm = TRUE)
    }
    
    forggplot_perc_1 <- do.call(rbind, forggplot_perc)
    forggplot_perc_1 <- data.frame(forggplot_perc_1, x = conc)
}

remove_peptides <- function(dataFrame, chans, reps, accessionID, chanNames, sequenceID, qualityID, 
    incPDofPD, PDofPD, removePeptides, modelType, incGeneFile, geneFile) {
    
    tempdat <- dataFrame
    
    
    # make matrix of descriptions and accessions to filter by common proteins
    
    accDesMat <- as.character(tempdat[as.character(tempdat[, accessionID]) == "", c(accessionID)])
    
    if (incPDofPD == TRUE) {
        data.merged <- tempdat[, c(chanNames, accessionID, sequenceID, qualityID, "pdofpd")]
        
    } else {
        data.merged <- tempdat[, c(chanNames, accessionID, sequenceID, qualityID)]
    }
    
    # need to make check if user has run the parameter function
    if (chanNames[1] != "Control_rep1") {
        channelnames <- standard_names(chans, reps, dataType = "intensity")
    }
    repindex <- rep(seq_len(reps), each = chans)
    totfal <- rep(FALSE, (chans + 3))
    
    if (incPDofPD == TRUE) {
        
        newdf <- cbind(data.merged[, seq_len(chans)], data.merged[, c(accessionID, sequenceID, qualityID, 
            PDofPD)])
        colnames(newdf) <- c(channelnames[repindex == 1], "Accession", "Sequence", "Quality", "Kd")
        newdf <- newdf[!is.na(rowSums(newdf[, seq_len(chans)])), ]
        newdf <- newdf[newdf$Quality <= 0.05, ]
        newdf <- data.frame(newdf, outliers = rep(0, length(newdf[, 1])), uniquePeps = rep(0, length(newdf[, 
            1])), addedVals = rep(0, length(newdf[, 1])), Kd = newdf$Kd)
        
        
    } else {
        newdf <- cbind(data.merged[, seq_len(chans)], data.merged[, c(accessionID, sequenceID, qualityID)])
        colnames(newdf) <- c(channelnames[repindex == 1], "Accession", "Sequence", "Quality")
        newdf <- newdf[!is.na(rowSums(newdf[, seq_len(chans)])), ]
        newdf <- newdf[newdf$Quality <= 0.05, ]
    }
    
    if (reps == 2) {
        channelnames <- paste("rep", rep(seq_len(reps), each = chans), "_", rep(paste("C", c("ontrol", 
            0:(chans - 2)), sep = ""), reps), sep = "")
        newdf2 <- cbind(data.merged[, (chans + 1):(2 * chans)], data.merged[, c(accessionID, sequenceID, 
            qualityID)])
        colnames(newdf2) <- c(channelnames[repindex == 2], "Accession", "Sequence", "Quality")
        newdf2 <- newdf2[!is.na(rowSums(newdf2[, seq_len(chans)])), ]
        newdf2 <- newdf2[newdf2$Quality <= 0.05, ]
        
        
        
        newdf <- data.frame(newdf, outliers = rep(0, length(newdf[, 1])), uniquePeps = rep(0, length(newdf[, 
            1])), addedVals = rep(0, length(newdf[, 1])))
        newdf2 <- data.frame(newdf2, outliers = rep(0, length(newdf2[, 1])), uniquePeps = rep(0, 
            length(newdf2[, 1])), addedVals = rep(0, length(newdf2[, 1])))
        
        
        common.proteins <- intersect(unique(newdf$Accession), unique(newdf2$Accession))
        
        newdf <- newdf[!is.na(match(newdf$Accession, common.proteins)), ]
        newdf2 <- newdf2[!is.na(match(newdf2$Accession, common.proteins)), ]
    } else {
        common.proteins <- unique(newdf$Accession)
    }
    
    if (removePeptides == FALSE) {
        
        if (reps == 1 & modelType == "sigmoid") {
            if (incPDofPD == TRUE) {
                prot1 <- unique(newdf$Accession)
                sumkd <- rep(0, length(prot1))
                protdf <- newdf[seq_len(length(prot1)), colnames(newdf) != "Sequence"]
                for (i in seq_len(length(prot1))) {
                  ## grep for total intensity () includes non unique peps
                  protdf[i, seq_len(chans)] <- apply(newdf[grep(prot1[i], newdf$Accession), seq_len(chans)], 
                    2, sum, na.rm = TRUE)
                  protdf$Kd[i] <- sum(newdf$Kd[grep(prot1[i], newdf$Accession)], na.rm = TRUE)
                  protdf$Accession[i] <- prot1[i]
                  
                  ## use == to get unique peptides per protein per repeat
                  protdf$uniquePeps[i] <- length(unique(unique(newdf$Sequence[newdf$Accession == 
                    prot1[i]])))
                  
                }
                
                protdf$uniquePeps[grep(";", protdf$Accession)] <- 0
                
                fc1 <- protdf[, 1]/protdf[, 2:chans]
                
                Kd <- protdf$Kd/protdf[, 1]
                
                fcprotdf <- data.frame(log2(fc1), Accession = protdf$Accession, uniquePepr1 = protdf$uniquePeps, 
                  uniquePepr2 = protdf$uniquePeps, Kd = Kd)
                fcprotdf <- fcprotdf[!is.na(rowSums(fcprotdf[seq_len((chans - 1))])), ]
                fcprotdf
                
            } else {
                
                prot1 <- unique(newdf$Accession)
                protdf <- newdf[seq_len(length(prot1)), colnames(newdf) != "Sequence"]
                
                for (i in seq_len(length(prot1))) {
                  ## grep for total intensity () includes non unique peps
                  protdf[i, seq_len(chans)] <- apply(newdf[grep(prot1[i], newdf$Accession), seq_len(chans)], 
                    2, sum, na.rm = TRUE)
                  protdf$Accession[i] <- prot1[i]
                  
                  ## use == to get unique peptides per protein per repeat
                  protdf$uniquePeps[i] <- length(unique(unique(newdf$Sequence[newdf$Accession == 
                    prot1[i]])))
                  
                }
                
                protdf$uniquePeps[grep(";", protdf$Accession)] <- 0
                
                fc1 <- protdf[, 1]/protdf[, 2:chans]
                
                
                fcprotdf <- data.frame(log2(fc1), Accession = protdf$Accession, uniquePepr1 = protdf$uniquePeps, 
                  uniquePepr2 = protdf$uniquePeps)
                fcprotdf <- fcprotdf[!is.na(rowSums(fcprotdf[seq_len((chans - 1))])), ]
                fcprotdf
                
            }
            
            
        } else {
            prot1 <- unique(newdf$Accession)
            protdf <- newdf[seq_len(length(prot1)), colnames(newdf) != "Sequence"]
            
            
            for (i in seq_len(length(prot1))) {
                ## grep for total intensity () includes non unique peps
                protdf[i, seq_len(chans)] <- apply(newdf[grep(prot1[i], newdf$Accession), seq_len(chans)], 
                  2, sum, na.rm = TRUE)
                
                protdf$Accession[i] <- prot1[i]
                
                ## use == to get unique peptides per protein per repeat
                protdf$uniquePeps[i] <- length(unique(c(unique(newdf$Sequence[newdf$Accession == 
                  prot1[i]]), unique(newdf2$Sequence[newdf2$Accession == prot1[i]]))))
                
            }
            
            protdf$uniquePeps[grep(";", protdf$Accession)] <- 0
            
            prot2 <- unique(newdf2$Accession)
            ### second data frame
            
            protdf2 <- newdf2[seq_len(length(prot2)), colnames(newdf2) != "Sequence"]
            
            for (i in seq_len(length(prot2))) {
                ## grep for total intensity () includes non unique peps
                protdf2[i, seq_len(chans)] <- matrixStats::colMedians(as.matrix(newdf2[grep(prot2[i], 
                  newdf2$Accession), seq_len(chans)]), na.rm = TRUE)
                protdf2$Accession[i] <- prot2[i]
                
                protdf2$uniquePeps[i] <- length(unique(newdf2$Sequence[newdf2$Accession == prot2[i]]))
            }
            
            protdf2$uniquePeps[grep(";", protdf2$Accession)] <- 0
            
            
            com.prot <- intersect(protdf$Accession, protdf2$Accession)
            
            fc1 <- protdf[, 1]/protdf[, 2:chans]
            fc2 <- protdf2[, 1]/protdf2[, 2:chans]
            
            
            fcprotdf <- data.frame(log2(fc1[match(com.prot, protdf$Accession), ]), log2(fc2[match(com.prot, 
                protdf2$Accession), ]), Accession = com.prot, uniquePepr1 = protdf$uniquePeps[match(com.prot, 
                protdf$Accession)], uniquePepr2 = protdf2$uniquePeps[match(com.prot, protdf2$Accession)])
            fcprotdf <- fcprotdf[!is.na(rowSums(fcprotdf[seq_len((reps * chans - 2))])), ]
            fcprotdf
        }
        
        
    } else {
        
        totpepdf <- NULL
        totpepdf2 <- NULL
        
        for (z in seq_len(length(common.proteins))) {
            temp <- newdf[newdf$Accession == common.proteins[z], ]
            temp2 <- newdf2[newdf2$Accession == common.proteins[z], ]
            # first step: check if all peptides are unique ...
            
            if (all.equal(grep(";", temp$Accession), integer(0)) == TRUE) {
                
                uniPeptides1 <- length(unique(temp$Sequence))
                
            } else {
                uniPeptides1 <- length(unique(temp$Sequence[-grep(";", temp$Accession)]))
            }
            
            if (all.equal(grep(";", temp2$Accession), integer(0)) == TRUE) {
                
                uniPeptides2 <- length(unique(temp2$Sequence))
                
            } else {
                
                uniPeptides2 <- length(unique(temp2[-grep(";", temp2$Accession)]))
            }
            
            ## add unique peptide column
            
            temp$uniquePeps <- uniPeptides1
            temp2$uniquePeps <- uniPeptides2
            
            
            tempPep <- intersect(unique(temp$Sequence), unique(temp2$Sequence))
            
            
            if (all.equal(tempPep, character(0)) == TRUE) {
                next
            }
            
            for (i in seq_len(length(tempPep))) {
                
                if (sum(temp$Sequence == tempPep[i]) != sum(temp2$Sequence == tempPep[i])) {
                  
                  dr <- peptide_match(temp[temp$Sequence == tempPep[i], ], temp2[temp2$Sequence == 
                    tempPep[i], ], chans)
                  dr1 <- dr$dr1
                  dr2 <- dr$dr2
                  
                } else {
                  
                  dr1 <- temp[temp$Sequence == tempPep[i], ]
                  dr2 <- temp2[temp2$Sequence == tempPep[i], ]
                  
                  
                }
                
                tempoindex <- rep(FALSE, nrow(dr1))
                for (j in seq_len(nrow(dr1))) {
                  percor <- stats::cor.test(log2(as.numeric(dr1[j, seq_len(chans)])), log2(as.numeric(dr2[j, 
                    seq_len(chans)])))
                  tempoindex[j] <- percor$estimate < 0.4
                }
                
                dr1[tempoindex, seq_len(chans)] <- NA
                dr2[tempoindex, seq_len(chans)] <- NA
                
                
                
                tempReplace <- dr1
                tempReplace2 <- dr2
                
                
                
                temp <- temp[temp$Sequence != tempPep[i], ]
                temp <- rbind(temp, tempReplace)
                
                temp2 <- temp2[temp2$Sequence != tempPep[i], ]
                temp2 <- rbind(temp2, tempReplace2)
                
            }
            
            temp <- temp[match(tempPep, temp$Sequence), ]
            temp2 <- temp2[match(tempPep, temp2$Sequence), ]
            
            totpepdf <- rbind(totpepdf, temp)
            
            
            totpepdf2 <- rbind(totpepdf2, temp2)
        }
        
        
        
        
        totpepdf <- totpepdf[!is.na(rowSums(totpepdf[, seq_len(chans)])), ]
        totpepdf2 <- totpepdf2[!is.na(rowSums(totpepdf2[, seq_len(chans)])), ]
        
        totpepdf <- totpepdf[totpepdf$addedVals == 0, ]
        totpepdf <- totpepdf[totpepdf$addedVals == 0, ]
        
        
        totpepdf$uniquePeps[grep(";", totpepdf$Accession)] <- 0
        
        
        pepframe <- data.frame(totpepdf[seq_len(length(common.proteins)), seq_len(chans)], totpepdf2[seq_len(length(common.proteins)), 
            seq_len(chans)], Accession = totpepdf$Accession[seq_len(length(common.proteins))], uniquePeps = totpepdf$uniquePeps[seq_len(length(common.proteins))])
        
        pepsum1 <- pepsum2 <- totpepdf[seq_len(length(common.proteins)), ]
        pepsum1 <- pepsum1[, -match(c("Sequence", "addedVals", "Quality", "outliers"), colnames(pepsum1))]
        pepsum2 <- pepsum2[, -match(c("Sequence", "addedVals", "Quality", "outliers"), colnames(pepsum2))]
        
        pepsum1$pepNum <- pepsum1$pepNum <- rep(0, length(common.proteins))
        
        
        colnames(pepsum2)[seq_len(chans)] <- channelnames[(chans + 1):(chans * reps)]
        
        for (i in seq_len(length(common.proteins))) {
            pepsum1[i, seq_len(chans)] <- colSums(totpepdf[grep(common.proteins[i], totpepdf$Accession), 
                seq_len(chans)])
            pepsum2[i, seq_len(chans)] <- colSums(totpepdf2[grep(common.proteins[i], totpepdf2$Accession), 
                seq_len(chans)])
            
            pepsum1$pepNum[i] <- nrow(totpepdf[grep(common.proteins[i], totpepdf$Accession), seq_len(chans)])
            pepsum2$pepNum[i] <- nrow(totpepdf2[grep(common.proteins[i], totpepdf2$Accession), seq_len(chans)])
            
            
            pepsum1$Accession[i] <- pepsum2$Accession[i] <- common.proteins[i]
            pepsum1$uniquePeps[i] <- totpepdf$uniquePeps[as.logical(match(totpepdf$Accession, common.proteins[i], 
                nomatch = FALSE))][1]
            
            
            
            
            
        }
        
        indexpepsum <- ((rowSums(pepsum1[seq_len(chans)]) != 0) + ((rowSums(pepsum2[seq_len(chans)])) != 
            0)) == 2
        pepsum1 <- pepsum1[indexpepsum, ]
        pepsum2 <- pepsum2[indexpepsum, ]
        
        
        
        fc1 <- pepsum1[, 1]/pepsum1[, 2:chans]
        fc2 <- pepsum2[, 1]/pepsum2[, 2:chans]
        
        
        fcprotdf <- data.frame(log2(fc1), log2(fc2), pepsum1$Accession, uniquePepr1 = pepsum1$uniquePeps, 
            uniquePepr2 = pepsum2$uniquePeps, num1 = pepsum1$pepNum, num2 = pepsum2$pepNum)
        
        
        fcprotdf
        
    }
    fcprotdf
}

normalize_data <- function(dataFrame, chans, reps, PD2, channelNames, incPDofPD, PDofPD, removePeptides, 
    dataType, modelType, incGeneFile = FALSE, conversionTable = NA, normaliseData = "median", accessionID = NA, 
    uniquePeptides = NA, organism = "H.sapiens") {
    
    inten <- dataFrame
    channelIndex <- seq_len((reps * chans))
    if (dataType != "intensity") {
        
        if (modelType == "sigmoid") {
            
            data_orig2 <- dataFrame
            
            ## if user has specified accession & Description
            if (PD2 == TRUE) {
                pattern <- "GN=(\\S+)"
                g_fromD1 <- stringr::str_extract(data_orig2$Description, pattern)
                gID_D1a <- stringr::str_split_fixed(g_fromD1, "GN=", n = 2)
                gID_D1a <- as.vector(gID_D1a[, 2])
                gID_D1 <- as.matrix(replace(gID_D1a, gID_D1a == "", "NA"))
                # Addition of the gene ID column
                data_orig2["geneID"] <- (gID_D1)
                
            } else {
                
                if (incGeneFile == FALSE) {
                  
                  tempacc <- data_orig2[, accessionID]
                  data_orig2 <- data_orig2[, (-accessionID)]
                  data_orig2$Accession <- tempacc
                  uniGene <- uniprotGene(organism)
                  uniGene$Gene.names[uniGene$Gene.names == ""] <- NA
                  GeneID <- uniGene$Gene.names[match(data_orig2$Accession, uniGene$Entry)]
                  GeneID <- make.names(GeneID, unique = TRUE)
                  
                  data_orig2$geneID <- GeneID
                } else {
                  
                  tempacc <- data_orig2[, accessionID]
                  data_orig2 <- data_orig2[, (-accessionID)]
                  data_orig2$Accession <- tempacc
                  uniGene <- conversionTable
                  GeneID <- uniGene$Gene.names[match(data_orig2$Accession, uniGene$Entry)]
                  GeneID <- make.names(GeneID, unique = TRUE)
                  data_orig2$geneID <- GeneID
                }
            }
            if (incPDofPD == TRUE) {
                if (PD2 == TRUE) {
                  data.merged <- data.frame(data_orig2[, channelNames], data_orig2$Accession, data_orig2$geneID, 
                    data_orig2$X..Unique.Peptides, data_orig2[, PDofPD])
                } else {
                  data.merged <- data.frame(data_orig2[, channelNames], data_orig2$Accession, data_orig2$geneID, 
                    data_orig2[, uniquePeptides], data_orig2[, PDofPD])
                }
                final.Names <- standard_names(chans, reps, dataType)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "Kd")
                tmp <- data.merged[, seq_len(chans)]
                tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
                countsNAs <- as.data.frame(apply(tf, 1, function(x) table(x)["TRUE"]))
                n_of_miss <- as.data.frame(as.numeric(stringr::str_replace_all(as.list(countsNAs[, 
                  1]), "NA", "0")))
                data.merged <- data.frame(data.merged, n_of_miss)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "Kd", 
                  "MissingVal")
                missing_val <- 0
                data.merged <- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
                # Specify the number of missing points.For zero missing point is ==0 filiter for 2 unique
                # peptides
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                if (normaliseData == "loess") {
                  data.merged <- data.frame((affy::normalize.loess(2^(data.merged[, channelIndex]))), 
                    Accession = data.merged$Accession, GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, 
                    Kd = 1/data.merged$Kd, MissingVal = data.merged$MissingVal)
                } else if (normaliseData == "median") {
                  data.merged <- data.frame(((2^(data.merged[, channelIndex]))/apply((2^(data.merged[, 
                    channelIndex])), 2, stats::median)), Accession = data.merged$Accession, GeneID = data.merged$GeneID, 
                    UniquePeps = data.merged$UniquePeps, Kd = 1/data.merged$Kd, MissingVal = data.merged$MissingVal)
                } else {
                  data.merged <- data.frame((2^(data.merged[, channelIndex])), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, Kd = 1/data.merged$Kd, 
                    MissingVal = data.merged$MissingVal)
                }
                data.merged
            } else {
                if (PD2 == TRUE) {
                  data.merged <- data.frame(data_orig2[, channelNames], data_orig2$Accession, data_orig2$geneID, 
                    data_orig2$X..Unique.Peptides)
                  
                } else {
                  data.merged <- data.frame(data_orig2[, channelNames], data_orig2$Accession, data_orig2$geneID, 
                    data_orig2[, uniquePeptides])
                }
                final.Names <- standard_names(chans, reps, "LFC")
                
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps")
                
                
                tmp <- data.merged[, seq_len(chans)]
                tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
                countsNAs <- as.data.frame(apply(tf, 1, function(x) table(x)["TRUE"]))
                n_of_miss <- as.data.frame(as.numeric(stringr::str_replace_all(as.list(countsNAs[, 
                  1]), "NA", "0")))
                data.merged <- data.frame(data.merged, n_of_miss)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "MissingVal")
                # 
                missing_val <- 0
                # 
                data.merged <- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
                # Specify the number of missing points.For zero missing point is ==0 filiter for 2 unique
                # peptides
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                if (normaliseData == "loess") {
                  
                  data.merged <- data.frame((affy::normalize.loess(2^(data.merged[, channelIndex]))), 
                    Accession = data.merged$Accession, GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, 
                    MissingVal = data.merged$MissingVal)
                  
                } else if (normaliseData == "median") {
                  
                  data.merged <- data.frame(((2^(data.merged[, channelIndex]))/apply((2^(data.merged[, 
                    channelIndex])), 2, stats::median)), Accession = data.merged$Accession, GeneID = data.merged$GeneID, 
                    UniquePeps = data.merged$UniquePeps, MissingVal = data.merged$MissingVal)
                  
                } else {
                  data.merged <- data.frame((2^(data.merged[, channelIndex])), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, MissingVal = data.merged$MissingVal)
                }
                data.merged
            }
            
        } else {
            
            tempdat <- dataFrame
            
            
            temp <- channelIndex
            nvec <- length(temp)
            
            data.merged <- tempdat[, channelNames]
            colnames(data.merged) <- standard_names(chans, reps, dataType)
            
            data.names <- c(colnames(data.merged), "Accession", "GeneID", "UniquePeps", "MissingVal")
            
            if (dataType == "FC") {
                data.merged[, temp] <- log2(data.merged[, temp])
            }
            
            if (PD2 == TRUE) {
                pattern <- "GN=(\\S+)"
                g_fromD1 <- stringr::str_extract(tempdat$Description, pattern)
                gID_D1a <- stringr::str_split_fixed(g_fromD1, "GN=", n = 2)
                gID_D1a <- as.vector(gID_D1a[, 2])
                gID_D1 <- as.matrix(replace(gID_D1a, gID_D1a == "", "NA"))
                
                Accession <- tempdat$Accession
                UniquePeps <- MissingVal <- tempdat[, grep("Unique", colnames(tempdat))]
                data.merged <- cbind(data.merged, Accession, GeneID = gID_D1, UniquePeps, MissingVal)
                
            } else {
                
                tempdat$Accession <- tempdat[, accessionID]
                if (incGeneFile == FALSE) {
                  
                  uniGene <- uniprotGene(organism)
                } else {
                  uniGene <- conversionTable
                  colnames(uniGene)[seq_len(2)] <- c("Entry", "Gene.names")
                }
                GeneID <- uniGene$Gene.names[match(tempdat$Accession, uniGene$Entry)]
                GeneID <- make.names(GeneID, unique = TRUE)
                
                Accession <- tempdat$Accession
                UniquePeps <- MissingVal <- tempdat[, uniquePeptides]
                
                data.merged <- cbind(data.merged, Accession, GeneID = GeneID, UniquePeps, MissingVal)
            }
            
            
            channelIndex <- seq_len((reps * chans))
            missing <- rowSums(is.na(data.merged[, channelIndex]))
            
            missing_val <- 0
            data.merged$MissingVal <- missing
            
            ## subset by missing
            data.merged <- data.merged[data.merged$MissingVal <= missing_val, ]
            data.merged <- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
            
            data.merged <- data.merged[data.merged$UniquePeps > 1, ]
            # 
            if (normaliseData == "loess") {
                
                data.merged <- data.frame(log2(affy::normalize.loess(2^(data.merged[, channelIndex]))), 
                  Accession = data.merged$Accession, GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, 
                  MissingVal = data.merged$MissingVal)
                
            } else if (normaliseData == "median") {
                
                data.merged <- data.frame(log2((2^(data.merged[, channelIndex]))/apply((2^(data.merged[, 
                  channelIndex])), 2, stats::median)), Accession = data.merged$Accession, GeneID = data.merged$GeneID, 
                  UniquePeps = data.merged$UniquePeps, MissingVal = data.merged$MissingVal)
                
            } else {
                data.merged <- data.frame(log2((2^(data.merged[, channelIndex]))), Accession = data.merged$Accession, 
                  GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, MissingVal = data.merged$MissingVal)
                
            }
            
            data.merged
            ## filter for 2 unique peptide
        }
        
    } else {
        
        ### intensities to protein done.. here no description so we need toi use accession and intermine
        
        tempdat <- inten
        
        if (incGeneFile == FALSE) {
            
            uniGene <- uniprotGene(organism)
        } else {
            uniGene <- conversionTable
            colnames(uniGene) <- c("Entry", "Gene.names")
        }
        
        if (removePeptides == TRUE) {
            
            Accession <- tempdat$pepsum1.Accession
        } else {
            Accession <- tempdat$Accession
        }
        uniGene$Gene.names[uniGene$Gene.names == ""] <- NA
        GeneID <- uniGene$Gene.names[match(Accession, uniGene$Entry)]
        GeneID <- make.names(GeneID, unique = TRUE)
        
        tempdat$GeneID <- GeneID
        
        UniquePeps <- tempdat$uniquePepr1
        
        if (reps == 1) {
            
            data.merged <- tempdat
            final.Names <- paste0("rep1_C", 0:(chans - 2))
            if (incPDofPD == TRUE) {
                
                data.merged <- data.frame(data.merged[, seq_len((chans - 1))], Accession = data.merged$Accession, 
                  GeneID = data.merged$GeneID, UniquePeps = data.merged$uniquePepr1, Kd = data.merged$Kd)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "Kd")
                
                tmp <- data.merged[, seq_len((chans - 1))]
                tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
                countsNAs <- as.data.frame(apply(tf, 1, function(x) table(x)["TRUE"]))
                n_of_miss <- as.data.frame(as.numeric(stringr::str_replace_all(as.list(countsNAs[, 
                  1]), "NA", "0")))
                data.merged <- data.frame(data.merged, n_of_miss)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "Kd", 
                  "MissingVal")
                # 
                missing_val <- 0
                # 
                data.merged <- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
                # Specify the number of missing points.For zero missing point is ==0 filter for 2 unique peptides
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                # 
                if (normaliseData == "loess") {
                  
                  data.merged <- data.frame((affy::normalize.loess(2^(data.merged[, seq_len((chans - 
                    1))]))), Accession = data.merged$Accession, GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, 
                    Kd = data.merged$Kd)
                } else if (normaliseData == "median") {
                  data.merged <- data.frame(((2^(data.merged[, seq_len((chans - 1))]))/apply((2^(data.merged[, 
                    seq_len((chans - 1))])), 2, stats::median)), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, Kd = data.merged$Kd)
                  
                } else {
                  data.merged <- data.frame((2^(data.merged[, seq_len((chans - 1))])), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps, Kd = data.merged$Kd)
                }
                
                
                data.merged
                
            } else {
                
                data.merged <- data.frame(data.merged[, seq_len((chans - 1))], Accession = data.merged$Accession, 
                  GeneID = data.merged$GeneID, UniquePeps = data.merged$uniquePepr1)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps")
                
                
                tmp <- data.merged[, seq_len((chans - 1))]
                tf <- as.data.frame(lapply(tmp, function(x) (is.na(x))))
                countsNAs <- as.data.frame(apply(tf, 1, function(x) table(x)["TRUE"]))
                n_of_miss <- as.data.frame(as.numeric(stringr::str_replace_all(as.list(countsNAs[, 
                  1]), "NA", "0")))
                data.merged <- data.frame(data.merged, n_of_miss)
                colnames(data.merged) <- c(final.Names, "Accession", "GeneID", "UniquePeps", "MissingVal")
                missing_val <- 0
                data.merged <- as.data.frame(data.merged[data.merged$MissingVal <= missing_val, ])
                # #Specify the number of missing points.For zero missing point is ==0 #filiter for 2 unique
                # peptides
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                if (normaliseData == "loess") {
                  
                  data.merged <- data.frame((affy::normalize.loess(2^(data.merged[, seq_len((chans - 
                    1))]))), Accession = data.merged$Accession, GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps)
                  
                } else if (normaliseData == "median") {
                  data.merged <- data.frame(((2^(data.merged[, seq_len((chans - 1))]))/apply((2^(data.merged[, 
                    seq_len((chans - 1))])), 2, stats::median)), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps)
                } else {
                  data.merged <- data.frame((2^(data.merged[, seq_len((chans - 1))])), Accession = data.merged$Accession, 
                    GeneID = data.merged$GeneID, UniquePeps = data.merged$UniquePeps)
                  
                }
                
                data.merged
            }
            
        } else {
            
            data.merged <- tempdat[, seq_len(((chans - 1) * reps))]
            
            data.merged <- tempdat[, seq_len(((chans - 1) * reps))]
            
            if (removePeptides == FALSE) {
                Accession <- tempdat$Accession
                MissingVal <- UniquePeps
                data.merged <- cbind(data.merged, Accession, GeneID = GeneID, UniquePeps, MissingVal)
                
                data.merged <- data.merged[!is.na(rowSums(data.merged[, seq_len(((chans - 1) * reps))])), 
                  ]
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                
                data.merged <- data.frame((log2(affy::normalize.loess(2^(data.merged[, seq_len(((chans - 
                  1) * reps))])))), Accession = data.merged$Accession, GeneID = data.merged$GeneID, 
                  UniquePeps = data.merged$UniquePeps)
                
                data.merged <- data.merged[!is.na(rowSums(data.merged[, seq_len(((chans - 1) * reps))])), 
                  ]
                data.merged
            } else {
                
                Accession <- tempdat$pepsum1.Accession
                data.merged <- cbind(data.merged, Accession, GeneID = GeneID, UniquePeps, num1 = tempdat$num1, 
                  num2 = tempdat$num2)
                
                
                data.merged <- data.merged[!is.na(rowSums(data.merged[, seq_len(((chans - 1) * reps))])), 
                  ]
                data.merged <- data.merged[data.merged$UniquePeps > 1, ]
                
                data.merged <- data.frame((log2(affy::normalize.loess(2^(data.merged[, seq_len(((chans - 
                  1) * reps))])))), Accession = data.merged$Accession, GeneID = data.merged$GeneID, 
                  UniquePeps = data.merged$UniquePeps, num1 = data.merged$num1, num2 = data.merged$num2)
                
                data.merged <- data.merged[data.merged$num1 > 1 | data.merged$num2 > 1, ]
                data.merged <- data.merged[!is.na(rowSums(data.merged[, seq_len(((chans - 1) * reps))])), 
                  ]
                data.merged
            }
            
        }
        
        
    }
    
    
    data.merged
    
}

fit_model <- function(dataFrame, chans, reps, dataType = "LFC", modelType = "sigmoid", sigmoidConc = NA, 
    PD2 = TRUE, incPDofPD = FALSE, PDofPD = NA) {
    
    if (modelType == "sigmoid") {
        if (dataType == "intensity") {
            nvec <- seq_len((chans - 1))
            nvec <- length(nvec)
            data.merged <- dataFrame
            
            conc <- sigmoidConc
            if (incPDofPD == TRUE) {
                
                final.Names <- paste0("rep1_C", 0:(chans - 2))
                pred.names <- paste0("predX", seq_len((chans - 1)))
                colnames(data.merged[, seq_len((chans - 1))]) <- final.Names
                data_merged_positives <- data.merged
                data_merged_positives2 <- ((1/data_merged_positives[, seq_len((chans - 1))])) * 100
                Reps_FC <- data.frame(data_merged_positives2, Accession = data_merged_positives$Accession, 
                  GeneID = data_merged_positives$GeneID, UniquePeps = data_merged_positives$UniquePeps, 
                  depletionConstant = data_merged_positives$Kd)
                
                
                ryegrass.m1 <- vector(mode = "list", length = nrow(Reps_FC))
                pvals <- vector(mode = "list", length = nrow(Reps_FC))
                stderr <- vector(mode = "list", length = nrow(Reps_FC))
                model_pred <- vector(mode = "list", length = nrow(Reps_FC))
                coeff_predicted <- vector(mode = "list", length = nrow(Reps_FC))
                for (i in seq_len(nrow(Reps_FC))) {
                  # maxIt and relTol to be user defined
                  ryegrass.m1[[i]] <- try(drc::drm(as.numeric(Reps_FC[i, seq_len((chans - 1))]) ~ 
                    as.numeric(conc), na.action = stats::na.omit, control = drc::drmc(constr = FALSE, 
                    errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06), fct = drc::LL.4(fixed = c(NA, 
                    NA, NA, NA), names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))), silent = TRUE)
                  
                }
                
                failed_sigm = 0
                for (i in seq_len(length(ryegrass.m1))) {
                  # checking_val if FALSE the model has failed to calculate the pval
                  checking_val <- try(is.numeric(vsn::coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]), 
                    silent = TRUE)
                  
                  if (checking_val == "TRUE") {
                    pvals[[i]] <- t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[13:16]))
                    colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                      "RB50Pval")
                    coeff_predicted[[i]] <- t(data.frame(vsn::coefficients(ryegrass.m1[[i]])))
                    colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                      "RB50Coef")
                    stderr[[i]] <- t(as.data.frame(summary(ryegrass.m1[[i]])$coefficients[5:8]))
                    colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                    
                    model_pred[[i]] <- predict(ryegrass.m1[[i]])
                  } else {
                    failed_sigm = failed_sigm + 1
                    
                    fit <- stats::lm(as.numeric(Reps_FC[i, seq_len((chans - 1))]) ~ poly(log10(conc), 
                      2))
                    # extract the pval
                    pval <- as.numeric(summary(fit)$coefficients[, 4])
                    pvals[[i]] <- t(as.data.frame(c(pval, "lm-fit:intercept.slope.quadratic")))
                    colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                      "RB50Pval")
                    stderr[[i]] <- data.frame(NA, NA, NA, NA)
                    colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                    coeff_predicted[[i]] <- data.frame(NA, NA, NA, NA)
                    colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                      "RB50Coef")
                    model_pred[[i]] <- as.numeric(stats::fitted(fit))
                  }
                }
                
                
                modelsReps <- data.frame(do.call(rbind.data.frame, lapply(model_pred, function(x) as.numeric(x))), 
                  Reps_FC$GeneID, do.call(rbind.data.frame, lapply(pvals, function(x) x)), do.call(rbind.data.frame, 
                    lapply(coeff_predicted, function(x) x)), do.call(rbind.data.frame, lapply(stderr, 
                    function(x) x)))
                colnames(modelsReps) <- c(pred.names, "GeneID", "SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                  "RB50Pval", "SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", "RB50Coef", "SlopeErr", 
                  "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                
                data_merged_2 <- merge.data.frame(modelsReps, Reps_FC, by = "GeneID")
                
                
                data_merged_2 <- data.frame(data_merged_2, Top_minus_min = data_merged_2$predX1 - 
                  data_merged_2[, paste("predX", (chans - 1), sep = "")])
                
                data_merged_2 <- data.frame(data_merged_2, correctedRB50 = (data_merged_2$RB50Coef * 
                  data_merged_2$depletionConstant))
                data_merged_2
                
            } else {
                
                final.Names <- paste0("rep1_C", 0:(chans - 2))
                pred.names <- paste0("predX", seq_len((chans - 1)))
                colnames(data.merged[, seq_len((chans - 1))]) <- final.Names
                data_merged_positives <- data.merged
                data_merged_positives2 <- ((1/data_merged_positives[, seq_len((chans - 1))])) * 100
                Reps_FC <- data.frame(data_merged_positives2, Accession = data_merged_positives$Accession, 
                  GeneID = data_merged_positives$GeneID, UniquePeps = data_merged_positives$UniquePeps)
                
                
                ryegrass.m1 <- vector(mode = "list", length = nrow(Reps_FC))
                pvals <- vector(mode = "list", length = nrow(Reps_FC))
                stderr <- vector(mode = "list", length = nrow(Reps_FC))
                model_pred <- vector(mode = "list", length = nrow(Reps_FC))
                coeff_predicted <- vector(mode = "list", length = nrow(Reps_FC))
                for (i in seq_len(nrow(Reps_FC))) {
                  # maxIt and relTol to be user defined
                  ryegrass.m1[[i]] <- try(drc::drm(as.numeric(Reps_FC[i, seq_len((chans - 1))]) ~ 
                    as.numeric(conc), na.action = stats::na.omit, control = drc::drmc(constr = FALSE, 
                    errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06), fct = drc::LL.4(fixed = c(NA, 
                    NA, NA, NA), names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))), silent = TRUE)
                  
                }
                
                failed_sigm = 0
                for (i in seq_len(length(ryegrass.m1))) {
                  # checking_val if FALSE the model has failed to calculate the pval
                  checking_val <- try(is.numeric(vsn::coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]), 
                    silent = TRUE)
                  
                  if (checking_val == "TRUE") {
                    pvals[[i]] <- t(as.data.frame(suppressWarnings(summary(ryegrass.m1[[i]])$coefficients[13:16])))
                    colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                      "RB50Pval")
                    coeff_predicted[[i]] <- t(data.frame(suppressWarnings(vsn::coefficients(ryegrass.m1[[i]]))))
                    colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                      "RB50Coef")
                    stderr[[i]] <- t(as.data.frame(suppressWarnings(summary(ryegrass.m1[[i]])$coefficients[5:8])))
                    colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                    
                    model_pred[[i]] <- suppressWarnings(predict(ryegrass.m1[[i]]))
                  } else {
                    failed_sigm = failed_sigm + 1
                    
                    fit <- stats::lm(as.numeric(Reps_FC[i, seq_len((chans - 1))]) ~ poly(log10(conc), 
                      2))
                    # extract the pval
                    pval <- as.numeric(summary(fit)$coefficients[, 4])
                    pvals[[i]] <- t(as.data.frame(c(pval, "lm-fit:intercept.slope.quadratic")))
                    colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                      "RB50Pval")
                    stderr[[i]] <- data.frame(NA, NA, NA, NA)
                    colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                    coeff_predicted[[i]] <- data.frame(NA, NA, NA, NA)
                    colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                      "RB50Coef")
                    model_pred[[i]] <- as.numeric(stats::fitted(fit))
                  }
                }
                
                
                modelsReps <- data.frame(do.call(rbind.data.frame, lapply(model_pred, function(x) as.numeric(x))), 
                  Reps_FC$GeneID, do.call(rbind.data.frame, lapply(pvals, function(x) x)), do.call(rbind.data.frame, 
                    lapply(coeff_predicted, function(x) x)), do.call(rbind.data.frame, lapply(stderr, 
                    function(x) x)))
                colnames(modelsReps) <- c(pred.names, "GeneID", "SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                  "RB50Pval", "SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", "RB50Coef", "SlopeErr", 
                  "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                
                data_merged_2 <- merge.data.frame(modelsReps, Reps_FC, by = "GeneID")
                
                data_merged_2 <- data.frame(data_merged_2, Top_minus_min = data_merged_2$predX1 - 
                  data_merged_2[, paste("predX", (chans - 1), sep = "")])
                
                data_merged_2
                
            }
            
            
        } else {
            
            nvec <- seq_len(chans)
            nvec <- length(nvec)
            data.merged <- dataFrame
            
            conc <- sigmoidConc
            final.Names <- standard_names(chans, reps, dataType = "LFC")
            pred.names <- paste0("predX", seq_len(chans))
            colnames(data.merged[, seq_len(chans)]) <- final.Names
            data_merged_positives <- data.merged
            data_merged_positives2 <- ((1/data_merged_positives[, seq_len(chans)])) * 100
            if (incPDofPD == TRUE) {
                Reps_FC <- data.frame(data_merged_positives2, Accession = data_merged_positives$Accession, 
                  GeneID = data_merged_positives$GeneID, UniquePeps = data_merged_positives$UniquePeps, 
                  depletionConstant = data_merged_positives$Kd, MissingVal = data_merged_positives$MissingVal)
            } else {
                Reps_FC <- data.frame(data_merged_positives2, Accession = data_merged_positives$Accession, 
                  GeneID = data_merged_positives$GeneID, UniquePeps = data_merged_positives$UniquePeps, 
                  MissingVal = data_merged_positives$MissingVal)
            }
            
            ryegrass.m1 <- vector(mode = "list", length = dim(data_merged_positives)[1])
            pvals <- vector(mode = "list", length = dim(data_merged_positives)[1])
            stderr <- vector(mode = "list", length = dim(data_merged_positives)[1])
            model_pred <- vector(mode = "list", length = dim(data_merged_positives)[1])
            coeff_predicted <- vector(mode = "list", length = dim(data_merged_positives)[1])
            for (i in seq_len(nrow(Reps_FC))) {
                # maxIt and relTol to be user defined
                ryegrass.m1[[i]] <- suppressWarnings(try(drc::drm(as.numeric(Reps_FC[i, seq_len(chans)]) ~ 
                  as.numeric(conc), na.action = stats::na.omit, control = drc::drmc(constr = FALSE, 
                  errorm = FALSE, noMessage = TRUE, maxIt = 1000, relTol = 1e-06), fct = drc::LL.4(fixed = c(NA, 
                  NA, NA, NA), names = c("Slope", "Lower Limit", "Upper Limit", "RB50"))), silent = TRUE))
                
            }
            failed_sigm = 0
            for (i in seq_len(length(ryegrass.m1))) {
                # checking_val if FALSE the model has failed to calculate the pval
                checking_val <- try(is.numeric(vsn::coefficients(ryegrass.m1[[i]])[["Slope:(Intercept)"]]), 
                  silent = TRUE)
                
                if (checking_val == "TRUE") {
                  pvals[[i]] <- t(as.data.frame(suppressWarnings(summary(ryegrass.m1[[i]])$coefficients[13:16])))
                  colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", "RB50Pval")
                  coeff_predicted[[i]] <- t(data.frame(vsn::coefficients(ryegrass.m1[[i]])))
                  colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                    "RB50Coef")
                  stderr[[i]] <- t(as.data.frame(suppressWarnings(summary(ryegrass.m1[[i]])$coefficients[5:8])))
                  colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                  
                  model_pred[[i]] <- suppressWarnings(predict(ryegrass.m1[[i]]))
                } else {
                  failed_sigm = failed_sigm + 1
                  
                  fit <- stats::lm(as.numeric(Reps_FC[i, seq_len(chans)]) ~ poly(log10(conc), 2))
                  # extract the pval
                  pval <- as.numeric(summary(fit)$coefficients[, 4])
                  pvals[[i]] <- t(as.data.frame(c(pval, "lm-fit:intercept.slope.quadratic")))
                  colnames(pvals[[i]]) <- c("SlopePval", "Lower_LimitPval", "Upper_LimitPval", "RB50Pval")
                  stderr[[i]] <- data.frame(NA, NA, NA, NA)
                  colnames(stderr[[i]]) <- c("SlopeErr", "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
                  coeff_predicted[[i]] <- data.frame(NA, NA, NA, NA)
                  colnames(coeff_predicted[[i]]) <- c("SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", 
                    "RB50Coef")
                  model_pred[[i]] <- as.numeric(stats::fitted(fit))
                }
            }
            
            modelsReps <- data.frame(do.call(rbind.data.frame, lapply(model_pred, function(x) as.numeric(x))), 
                Reps_FC$GeneID, do.call(rbind.data.frame, lapply(pvals, function(x) x)), do.call(rbind.data.frame, 
                  lapply(coeff_predicted, function(x) x)), do.call(rbind.data.frame, lapply(stderr, 
                  function(x) x)))
            colnames(modelsReps) <- c(pred.names, "GeneID", "SlopePval", "Lower_LimitPval", "Upper_LimitPval", 
                "RB50Pval", "SlopeCoef", "Lower_LimitCoef", "Upper_LimitCoef", "RB50Coef", "SlopeErr", 
                "Lower_LimitErr", "Upper_LimitErr", "RB50Err")
            data_merged_2 <- merge.data.frame(modelsReps, Reps_FC, by = "GeneID")
            
            data_merged_2 <- data.frame(data_merged_2, Top_minus_min = data_merged_2$predX1 - data_merged_2[, 
                paste("predX", chans, sep = "")])
            
            if (incPDofPD == TRUE) {
                
                data_merged_2 <- data.frame(data_merged_2, correctedRB50 = (data_merged_2$RB50Coef * 
                  data_merged_2$depletionConstant))
            } else {
                data_merged_2 <- data.frame(data_merged_2)
            }
            data_merged_2
            
        }
        
    } else {
        if (dataType == "intensity") {
            nvec <- seq_len(((chans - 1) * (reps - 1)))
            nvec <- length(nvec)
            data.merged <- dataFrame
            conc <- rep(0:(chans - 2), times = reps)
            
        } else {
            nvec <- seq_len(((chans) * (reps)))
            nvec <- length(nvec)
            data.merged <- dataFrame
            conc <- rep(0:(chans - 1), times = reps)[seq_len((chans * reps))]
            
        }
        
        design <- stats::model.matrix(~poly(conc, 2))
        colnames(design) <- c("Intercept", "Slope", "Quadratic")
        
        fit <- limma::lmFit(data.merged[, seq_len(length(conc))], method = "ls", design = design)
        fit <- limma::eBayes(fit)
        res <- limma::topTable(fit, coef = "Slope", number = nrow(data.merged), adjust = "BH")  #pval for the slope
        res2 <- limma::topTable(fit, coef = 1, number = nrow(data.merged), adjust = "BH")  #pval for the intercept
        res3 <- limma::topTable(fit, coef = "Quadratic", number = nrow(data.merged), adjust = "BH")  #pval for the quadratic term ()
        
        # add the pvalues to the dataframe
        tmp_1 <- cbind(data.merged[rownames(res), ], res)
        
        tmp_2 <- cbind(data.merged[rownames(res2), ], res2)
        
        tmp_3 <- cbind(data.merged[rownames(res3), ], res3)
        
        tobeselected <- merge.data.frame(tmp_1, tmp_2, by = "Accession")
        tobeselected <- merge.data.frame(tobeselected, tmp_3, by = "Accession")
        
        selectnames <- c(paste0(colnames(data.merged)[seq_len(length(conc))], ".x"), "logFC.x", "AveExpr.x", 
            "P.Value", "adj.P.Val", "P.Value.x", "adj.P.Val.x", "P.Value.y", "adj.P.Val.y", "Accession", 
            "GeneID.x", "UniquePeps")
        # can create a replacement names function here ....
        data.merged <- tobeselected[, match(selectnames, colnames(tobeselected))]
        nam <- standard_names(chans, reps, dataType)
        if (dataType == "intensity") {
            nam <- nam[-seq(1, reps * chans, by = chans)]
        }
        colnames(data.merged)[seq_len(length(nam))] <- nam
        colnames(data.merged) <- c(nam, "logFC", "AveExpr", "P.Value_slope", "adj.P.Val_slope", "P.Value_intercept", 
            "adj.P.Val_intercept", "P.Value_quadratic", "adj.P.Val_quadratic", "Accession", "GeneID", 
            "UniquePeps")
        data.merged
    }
}
