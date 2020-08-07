#' a convient wrapper for creating plot for pedigree data
#' 
#' @param ped a pedigree that contains \code{ID, Sex, MotherID, FatherID, isProband, CurAge}
#' @param annot.cancers the cancers shortnames to display. When set to default /code{NULL}, the top four cancers will be displayed.
#' @param annot.feature ONE feature that we would get annotation from pedigree, one of the choices from \code{c("ethic","Twins","CurAge","race")}
#' @export
plotPed <- function(ped, annot.cancers="all", annot.features=NULL, cex=0.8, branch=1, ...) {
  
  if (FALSE) {
    annot.cancers <- c("BC","OC")
  }
  
  # change default to annotate all cancer ages,
  # when aff_age is missing, use "NA" to represent
  
  base_fam <- with(ped, {
    # in kinship2, 1=male, 2=female, 3=unknown
    # in our pedigree, Sex=ismale
    
    Sex <- ifelse( Sex == 0, 2, Sex)
    Sex <- ifelse(is.na(Sex), 3, Sex)
    
    # change NA code
    MotherID = ifelse(MotherID == -999, NA, MotherID)
    FatherID = ifelse(FatherID == -999, NA, FatherID)
    
    cbind(ID, FatherID, MotherID, Sex, isDead)
  })
  colnames(base_fam) <- c("id", "dadid", "momid", "sex", "censor")
  fam_vec <- kinship2::makefamid(id = base_fam[,"id"], 
                                 father.id = base_fam[,"dadid"],
                                 mother.id = base_fam[,"momid"])
  
  # if contains diconnected families
  if (length(unique(fam_vec)) > 1) {
    rlang::abort("Pedigree contains disconnected families, please first double-check by kinship2::makefamid", level = "DisconnectFamily")
  } 
  
  # if we have added missing parents
  if (!is.null(ped$add) && !all(ped$add == 0)) {
    rlang::warn("Missing parents are included in the plot", level = "PlotAddedParents")
    add_id <- ped$ID[ped$add == 1]
  }
  else {
    add_id <- NULL
  }
  
  cancers <- PanelPRO:::.getCancersFromFam(ped)
  if (annot.cancers != "all") {
    annot.cancers <- intersect(cancers, annot.cancers)
    affm <- as.matrix(ped[, paste0("isAff", annot.cancers)])
  }
  else {
    annot.cancers <- cancers
    affm <- as.matrix(ped[, paste0("isAff", cancers)])
  }
  
  affm[affm == -999] <- NA
  #affm[is.na(affm)] <- -1
  
  
  id2 <- NULL
  if (!is.null(annot.cancers)) {
    annot.cancers <- intersect(paste0("Age", PanelPRO:::.getCancersFromFam(ped)), colnames(ped))
    annot.cancers <- gsub("Age","",annot.cancers)
    Agecancers <- ped[, paste0("Age", annot.cancers)]
    agelist <-  list()
    for (ac in annot.cancers){
      ori <- ped[[paste0("Age", ac)]] 
      aff <- ped[[paste0("isAff", ac)]] 
      ori <- ifelse(aff, ori, NA)
      ori[is.na(ori) | ori == -999] <- ""
      ori <- ifelse(ori != "", paste(ac, ori, sep = ":"), "")
      agelist[[ac]] <- ori
    }
    id2 <- do.call(paste, c(agelist, list(sep = "\n")))
    #id2 <- ifelse(id2 == ",", "", id2 )
  }
  
  pedAll <- kinship2::pedigree(id = base_fam[,"id"],
                               dadid = base_fam[,"dadid"],
                               momid = base_fam[,"momid"],
                               sex = base_fam[,"sex"],
                               affected = affm,
                               status = base_fam[,"censor"])
  
  
  if (!is.null(annot.features)) {
    annot.features <- intersect(c("Twins", "ethic", "CurAge", "race"), annot.features)
    featuresdf <- ped[annot.features]
    idf <- lapply(1:nrow(ped), function(i){
      mark <- substr(featuresdf[i,], 1, 3)
      annot.features <- toupper(substr(annot.features, 1,1))[!is.na(mark)]
      mark <- mark[!is.na(mark)]
      paste(annot.features, mark, sep = ":")
    })
    
    id2 <- paste0(id2, sapply(idf, paste0, collapse = "\n"))
  }
  probands  = which(ped$isProband == 1)
  
  # infer the size
  plist <- kinship2::align.pedigree(pedAll)
  sinf <- dim(plist$nid)
  # dev.new(width = 0.2*sinf[2], height = 2*sinf[1])
  plot_backend(pedAll, annot = id2, 
               feature.name = annot.features, 
               which.proband = probands, add.ids= add_id, ...)
  # on.exit(dev.off())
  invisible(pedAll)
  
}