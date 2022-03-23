#' Bam files preparation for EventPointer
#'
#' Prepares the information contained in .bam files to be analyzed by EventPointer
#'
#' @param Samples Name of the .bam files to be analyzed (Sample1.bam,Sample2.bam,...,etc).
#' @param SamplePath Path where the bam files are stored.
#' @param Ref_Transc Reference transcriptome used to name the genes found in bam files. Options are: Ensembl, UCSC or GTF.
#' @param fileTransc Path to the GTF reference transcriptome ff Ref_Transc is GTF.
#' @param cores Number of cores used for parallel processing.
#' @param Alpha Internal SGSeq parameter to include or exclude regions
#'
#' @return SGFeaturesCounts object. It contains a GRanges object with the corresponding elements to build
#' the different splicing graphs found and the counts related to each of the elements.
#'
#' @examples
#' \dontrun{
#'  # Obtain the samples and directory for .bam files
#'
#'    BamInfo<-si
#'    Samples<-BamInfo[,2]
#'    PathToSamples <- system.file('extdata/bams', package = 'SGSeq')
#'    PathToGTF<-paste(system.file('extdata',package='EventPointer'),'/FBXO31.gtf',sep='')
#'
#'   # Run PrepareBam function
#'    SG_RNASeq<-PrepareBam_EP(Samples=Samples,
#'                             SamplePath=PathToSamples,
#'                             Ref_Transc='GTF',
#'                             fileTransc=PathToGTF,
#'                             cores=1)
#' }
#' @export


PrepareBam_EP <- function(Samples, 
    Ref_Transc = "Ensembl", fileTransc = NULL, 
    cores = 1, which=NULL, min_junction_count = 5, max_complexity = 100) {
    # Event Pointer for RNASeq Data
    cat("Preparing BAM files for EventPointer...")
    
    Bam_Info <- Samples

    
    cat("\n Obtaining Reference Transcriptome...")
    
    stopifnot(Ref_Transc == "Ensembl" | Ref_Transc == 
        "UCSC" | Ref_Transc == "GTF")
    
    if (Ref_Transc == "Ensembl") {
        TxDb <- makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
            dataset = "hsapiens_gene_ensembl", 
            host = "grch37.ensembl.org")
    } else if (Ref_Transc == "UCSC") {
        TxDb <- makeTxDbFromUCSC(genome = "hg19", 
            tablename = "knownGene")
    } else if (Ref_Transc == "GTF") {
        stopifnot(!is.null(fileTransc))
        
        TxDb <- GenomicFeatures:::makeTxDbFromGFF(file = fileTransc, 
            format = "gtf", dataSource = "External Transcriptome")
    } else {
        stop("Unknown Reference Transcriptome")
    }
    
    
    # Steps for the Reference Transcriptome
    
    # Convert the TxDb to Features
    # (GenomicFeatures)
    TxF_Ref <- convertToTxFeatures(TxDb)
    
    # Steps for bam files
    
    cat("Done")
    
    cat("\n Predicting Features from BAMs...")
    
    # Predict TxFeatures from the input bam
    # files
    TxF_mod <- SGSeqMod:::predictTxFeatures(Bam_Info, cores = cores, min_junction_count = min_junction_count, max_complexity = max_complexity)
    # TxF_mod <- SGSeqMod:::predictTxFeatures(Bam_Info, cores = 4, min_junction_count = 0, max_complexity = 100)
    
    # save(TxF_mod, file="prediction_reales_cx.RData")
    # load("prediction_reales_cx.RData")
    # TxF_mod <- SGSeqMod:::predictTxFeatures(Bam_Info[29,], which = which[which(seqnames(which) == "GL000220.1")][1] , cores = 1, min_junction_count = 0, max_complexity = 100)
    
    # Convert predicted Features to Splicing
    # Graph
    SgF_mod <- SGSeqMod:::convertToSGFeatures(TxF_mod)
    
    # Get the reads for each subexon and
    # junction

    SgFC_Mod <- SGSeqMod:::getSGFeatureCounts(Bam_Info, 
                                          SgF_mod, cores = cores)
    # save(SgFC_Mod, file="counts_reales_cx.RData")
    # load(file="counts_reales_cx.RData")
    # Relate the SG Features with the
    # Reference Transcriptome
    seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(SgFC_Mod)
    SgFC_Mod <- annotate(SgFC_Mod, TxF_Ref)


    return(SgFC_Mod)
}
