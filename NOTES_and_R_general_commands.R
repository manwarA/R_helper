library(stringr)
library(dplyr)

# This is for general history and commands list
#===================================
# General,  
#===================================
# Check doParallel package

# data.table fread is faster than regular read.csv
file <- data.table::fread("file.csv",
						  sep = "\t", 
						  blank.lines.skip=TRUE, 
						  header = TRUE)

# Run this command in case many background files are open e.g. plots or pdf files
dev.off()

# select max value between two columns; pmax is parallel maximum
common1$log10pvalue <- pmax(common1$log10pvalue.x, common1$log10pvalue.y)

# to chek the type of individual coloumn types
str(data_frame)

# Evat data in survival analysis should be numeric, lol. The factor will create problem,
# and to convert factor into numeric, first convert into character then into numeric, like -- as.numeric(as.character((data)). 
# Sometimes, conversion is not this straightforward.

# to check the dim of multiple dfs in a list.
lapply(result, dim)

# version of R-base and other relevant information
version

# check where lib are installed
.libPaths()

# check the common element between two vector
intersect(Up, Dn)
# "Metazoa_SRP" "Y_RNA" 

# reading excel files
# library needed (readxl)
readxl::read_excel

# If the library is installed, you can access any function from that lib using double colons (::), for example:
# readxl::read_excel

# list function in R package
ls("package:enrichR")

# In case, log2fc values of 2 diff gene are same, small increase small increment
#it will add 0.00001 to the second gene of the same log2fc value
mydata <- mydata %>%
		group_by(log2fc) %>% 
		mutate(log2fc2 = log2fc + seq(0, by=0.00001, length.out= n()))

# convert "character to numeric" in data frame
raw.ad[] <- sapply(raw.ad, as.numeric) ; it preserves both col names and row names, plus it requires [] on left side

# drop columns that have 1 factor level
df[sapply(df, nlevels) > 1]

# compare all files in the list
outer(allFiles, allFiles, Vectorize(all.equal))

#===================================
# multiple word replcement
#===================================
library(stringr)
text <- c("just writing an illustrative text for example,",
          "but it has some different text in each sentence,",
          "so I am just gonna replicate it in here.")

wrong_words <-  c("illustrative text", "it in here", "for example")
new_words <- c("illustrative_text", "it_in_here", "for_example")
for (i in seq_along(wrong_words)) {
    text <- str_replace_all(text, wrong_words[i], new_words[i])
}
text

#===================================
# No idea
#===================================
# string manipulation, regrex,  what for?
samples <- gsub("(?<=QC1)[^;]*", "", samples, perl = TRUE)
samples2 <- grepl("Unshared", samples2)
sum(grepl("LateStageTumor_pool.QC1", samples))
sum(grepl("Unshared.", samples))
samples2 <- grep("Unshared.", samples, value =T)
samples2 <- str_remove(samples2, "Unshared.Log.QC2_LateStageTumor_pool.QC1")
samples2 <- samples2[nzchar(samples2)] # nzchar is a fast way to find out if elements of a character vector are non-empty strings.


# string split or you can remove the remaining part from ENSEMBL name, the version of ensembl id such as ENSG0000000000012.3. Remove .3, otherwise megeing will be difficult.
rownames(df) <- sub("\\..*", "", rownames(df)) # gsub()

# strsplit use case
# If you need to extract the first (or nth) entry from each split, use:
word <- c('apple-orange-strawberry','chocolate')
sapply(strsplit(word,"-"), `[`, 1)
#[1] "apple"     "chocolate"

# Or faster and more explictly:
vapply(strsplit(word,"-"), `[`, 1, FUN.VALUE=character(1))
#[1] "apple"     "chocolate"

# Both bits of code will cope well with selecting whichever value in the split list, and will deal with cases that are outside the range:
vapply(strsplit(word,"-"), `[`, 2, FUN.VALUE=character(1))
#[1] "orange" NA  

#===================================
# t-test rowwise (Welch t-test) and STATISTICAL tests
#===================================
# by default, R performs Welch Two Sample t-test
x$stat <- sapply(1:nrow(x), function(i) 
	t.test(	as.numeric(as.character(unlist(x[i,2:4]))), 
			as.numeric(as.character(unlist(x[i,8:10])))
			)[c("p.value")])

# apply Wilcoxon test (this does not assume normal distribution)
wilcox.test(as.numeric(common2[1, 2:193]), random_dist, paired = F, alternative = "two.sided")

#===================================
# Limma for gene expression
#===================================
# limma for microarray data analysis
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# limma is for continuous data while the DESeq2 and EdgeR are for count data
# Moreover, in microarray, the Robust Muti-array analysis (RMA) converts the intensities into log2 form, making it easy to find log2FC
# and simply running the t test, the significance can be calculated

# calculate normalizing factor
d0 <- calcNormFactors(d0)

# Filter low-expressed genes
cutoff <- 1 # can be vary
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
# number of genes left
dim(d) 

plotMDS(d, col = as.numeric(group)) # col is the type of group you want to check the contrast/group

# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
mm <- model.matrix(~0 + group)

# OR create design matrix				 
design_limm <- model.matrix(~ factor(mapping3$sampleType))
# The above specifies a model where each coefficient corresponds to a group mean

# plot mean-variance trend
y <- voom(d, mm, plot = T)

fit <- lmFit(y, design_limm OR mm)
ebayes <- eBayes(fit)

lod <- -log10(ebayes[["p.value"]][, 2])
mtstat <- ebayes[["t"]][, 2]

# top significant genes
tab <- topTable(ebayes, coef=2, adjust="fdr", n=10)

# in case, new contrast has to be analyzed, just create new contrast and put that one in formula, rest is same
contr <- makeContrasts(groupI5.6 - groupC.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)				 

#===================================
# GEO related
#===================================						 
# create expression set for GEOquery 
# it was a test, but the GEO data is very unpredictable/messsy
# making it really hard to convert it into ExpressionSet

# for annotating hgu133plus2 Affymatrix data set
# BiocManager::install("hgu133plus2probe") for annotation
# BiocManager::install("gcrma")
# BiocManager::install("hgu133plus2cdf")
# BiocManager::install("hgu133plus2.db")

# for geo data
gset <- getGEO("GSE12056", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# mapping between human genome and illumine platforms (reference: https://github.com/AlexsLemonade/refinebio/issues/232)
# I've just gone through the human platforms. Here's a Bioconductor package to platform name mapping that I think should work:
illuminaHumanv1.db: Illumina Human-6 v1.0, Illumina HumanRef-8 v1.0
illuminaHumanv2.db: Illumina Human-6 V2.0, Illumina HumanRef-8 v2.0
illuminaHumanv3.db: Illumina HumanHT-12 V3.0, Illumina HumanRef-8 v3.0, Illumina HumanWG-6 V3.0
illuminaHumanv4.db: Illumina HumanHT-12 V4.0

#get data from GEO
geo_data <- getGEO(GEO = "GSE12456",
				   destdir = dir_path,
				  GSEMatrix = T, 
				  getGPL = F) # getGPL can be heavy
# read local data, series_matrix file
series_mat <- getGEO(filename = file_path,
					GSEMatrix = True)
				 
eset <- ExpressionSet(assayData = as.matrix(mat),
                      phenoData =  Biobase::AnnotatedDataFrame(phenoData(df[[1]])))

phenoData <- new("AnnotatedDataFrame",
                 data=pData(test2[[1]])) 
                 #varMetadata=metadata)                      

exampleSet <- ExpressionSet(assayData=as.matrix(mat),
                               phenoData=phenoData,
                               experimentData=experimentData(test2[[1]]),
                               annotation="hgu95av2")

#=========================================
# Methylation data analysis; manually downloaded files list
#=========================================

# all 146 files that were downloaded as the Normal RNA seq data also 
# contains the matched tumor RNA seq data
# convert meta data into factor for linear modelling, 

# minfi package was used to preprocess the data and get the beta values
library(minfi)

rgSet <- read.metharray.exp("TCGA_data/methylation_green_red/00_dataCombine/")

phenoData <- DataFrame("samples" = sampleNames(rgSet))
phenoData$group <- ifelse(grepl("normal", phenoData$samples), "Normal", "Tumor")
rownames(phenoData) <- phenoData$samples
phenoData <- phenoData[sampleNames(rgSet), ]
match(phenoData$samples, sampleNames(rgSet))
head(phenoData)

# adding phenotype data to methylation data (rgSet)
pData(rgSet) <- phenoData

# (1) quality control; needs to done using more than one matric

Mset <- preprocessRaw(rgSet) # nothing done, but the output is MethylSet
qc_bumphunter <- getQC(Mset) # class 'MethylSet' or 'GenomicMethylSet' required
plotQC(qc_bumphunter)

# (2) quality control using detection p-values; p-value > 0.05 should be avoided
detP <- detectionP(rgSet)

# examine mean detection p-values across all samples using bar-plot to identify any failed samples
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")

# (3) quality Control: The overall density distribution of Beta values for each sample is another useful
# metric to determine sample quality. Usually, one would expect to see most Beta 
# values to be either close to 0 or 1, indicating most of the CpG sites in the sample 
# are unmethylated or methylated. The densityPlot function plots these distribution 
# for each sample.

phenoData <- pData(MSet)
densityPlot(Mset, sampGroups = phenoData$group)

# Normalization should be done using Funnorm method if the study comprises of cancer, normal samples.
grSet <- preprocessFunnorm(rgSet, ratioConvert = T)

# Compare with the unnormalized data to visualize the effect of the normalization. 
# First a comparison of the Beta distributions for the different probe designs. 
# This will give an indication of the effectiveness of the within-array normalization.

par(mfrow=c(1,1))
# Plot distributions prior to normalization for sample 1
plotBetasByType(Mset[,1], main="Raw")

# The normalized object is a GenomicRatioSet which does not contain
# the necessary probe info, we need to extract this from the MethylSet first.
typeI <- getProbeInfo(Mset, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(Mset, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))

# Now plot the distributions of the normalized data for sample 1
plotBetasByType(getBeta(grSet)[,1], probeTypes = probeTypes, main="Normalized",)

# Does it look like the normalization brought the distributions closer to each other? 
# Now let’s see how the between-array normalization worked.
# visualise what the data looks like before and after normalization

library("RColorBrewer")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=phenoData$group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(phenoData$group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(grSet), sampGroups=phenoData$group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(phenoData$group)),
       text.col=brewer.pal(8,"Dark2"))


# Filtering
# Poor performing probes can obscure the biological signals in the data and are 
# generally filtered out prior to differential methylation analysis. 
# As the signal from these probes is unreliable, by removing them we perform
# fewer statistical tests and thus lower the multiple testing penalty. before filtering,
# ensure that probes are in the same order in the mSetSq and detP objects
detP <- detectionP(rgSet)
detP <- detP[match(featureNames(grSet), rownames(detP)),]

# remove any probes that have failed in one or more samples; this next line
# checks for each row of detP whether the number of values < 0.01 is equal
# to the number of samples (TRUE) or not (FALSE)
keep <- rowSums(detP < 0.01) == ncol(grSet)
table(keep)
# Subset the GenomicRatioSet
grSetFlt <- grSet[keep,]
grSetFlt

# Because the presence of short nucleotide polymorphisms (or SNPs) inside the 
# probe body or at the nucleotide extension can have important consequences on 
# the downstream analysis, minfi offers the possibility to remove such probes.

grSetFlt <- dropLociWithSnps(grSetFlt)
grSetFlt

# Once the data has been filtered and normalised, it is often useful to re-examine 
# the MDS plots to see if the relationship between the samples has changed. 
# From the new MDS plots it is apparent that much of the inter-individual variation 
# has been removed as this is no longer the first principal component, likely due to 
# the removal of the SNP-affected CpG probes. However, the samples do still cluster 
# by individual in the second dimension and thus a factor for individual should 
# still be included in the model.

library('limma')
library("Gviz")
library("DMRcate")
library("DMRcatedata")
library("stringr")
library("mCSEA")

# set up the grouping variables and colours
pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(phenoData$group))]
names(groups) <- levels(factor(phenoData$group))
cols <- groups[as.character(factor(phenoData$group))]

par(mfrow=c(1,2))
plotMDS(getM(grSetFlt), top=1000, gene.selection="common",
        col=pal[factor(phenoData$group)], cex=0.8)
legend("right", legend=levels(factor(phenoData$group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(grSetFlt), top=1000, gene.selection="common",
        col=pal[factor(phenoData$group)])
legend("right", legend=levels(factor(phenoData$group)), text.col=pal,
       cex=0.7, bg="white")


# create cluster of probes to be analyzed together
# grSet_with_beta <- ratioConvert(grSet)
# grSet_with_beta_grange <- mapToGenome(grSet_with_beta)
# beta.values <- getBeta(grSet_with_beta)
# dim(beta.values) # 48551 49

phenoData$group <- relevel(as.factor(phenoData$group), ref = "Normal")

#
cluster_bumphunter <- clusterMaker(grSet, maxGap = 300)

# using Bumphunter
design_for_bumphunter <- model.matrix(~0 + group, data = phenoData)

parallel_cluster <- makeCluster(6)
bumphunter_output <- bumphunter(grSet, 
								design = design_for_bumphunter,
                                coef = 2, 
								pickCutoff = TRUE, 
								nullMethod=c("permutation"), 
                                B = 49, 
								type = "Beta")
stopCluster(parallel_cluster)

#=========================================
# SeSaME package to analyze differential methylation
#=========================================
library(sesame)
library(BiocParallel)

# need to add the command sequence, however, it is similar to any model fitting method;
# (1) create model matrix
# (2) linear model fitting
# (3) creation of contrasts (optional)
# (4) analysis
				 
#===================================
# Process multiple files
#===================================
# Make a function to process each file, file name or identifier has to be appended to the respective columns. 
# after that, using sapply, this function can be used to all the files.

allFiles <- lapply(listOfFiles, function(x) readr::read_tsv(x,
                                                col_names = T,
                                                skip_empty_rows = T,
                                                trim_ws = TRUE)

# pattern based matching and retreiving the data
output1 <- Sys.glob("NCC_*_Proteome_KU_*\\OUTPUT\\")
paths <- list.files(output1, 
					pattern= glob2rx("*summed_tum_normal_refine*.csv$*"),
                    full.names=T, 
					recursive=T)
				 
processFile <- function(f) {
  bname = strsplit(basename(f), '_')[[1]][1]
  df = data.table::fread(f, 
						 sep = "\t", 
						 blank.lines.skip=TRUE, 
						 header = TRUE, 
                         select = c(1,6:9), 
						 data.table = TRUE)
	# create and assign new col names for each file. 
  colnames(df) = c('uniprot', 
				   paste0(bname, '_1'), 
				   paste0(bname, '_2'), 
                   paste0(bname, '_3'), 
				   paste0(bname, '_4'))
  df_list = append(df_list, list(df))
  df_list 
                        }

# Apply the function to all files.
result <- sapply(paths, processFile)

# another way to read multiple files
reader.maf.files <- function(fpath, numb){
    bname = strsplit(basename(fpath), split = ".", fixed = TRUE,)[[1]][1]
    file = read.table(fpath, sep = "\t", header = TRUE)
    file = file[col.to.keep]
    colnames(file)[1] <- bname
    colnames(file)[3:length(file)] = paste0(colnames(file)[3:length(file)], "_", counter)
    counter <<- counter + 1  # "<<-" to update values outside of function.
    return(file)
}

counter = 1
maf.luad <- lapply(files, reader.maf.files)
				 
#===================================
# Merging
#===================================
# Multimerge, merge multiple dfs, if the number of dfs are large and the size is also big, it can run out of memory.
# allow.cartesian = TRUE is default in DF, while in data.table merge, it is FALSE, 
# you need to set it explicitly to proceed merging of data.tables.

# base R solution
merged_output <- Reduce(function(x, y) merge(x, y, all.x = TRUE, by = c("uniprot"), 
                                            allow.cartesian=TRUE), result[1:40])

# purr based solution
list_of_data %>% purrr::reduce(left_join, by = "row_names") # for purr based solution to merge multiple dataframes; but it needs a unique col name (e.g. row_names) in each df.

# data table approach,let see how efficient it is
library(data.table)
setkey(result_merge2, "uniprot")
df2_3 <- as.data.table(merge2)[as.data.table(merge1), on = "uniprot"]# allow.cartesian=TRUE ]

#result_merge <- result_merge[rowSums(is.na(result_merge[, 2:ncol(result_merge)])) == 0, ]
#result_merge <- result_merge %>%  dplyr::select(-starts_with("gene"))

## NOTE: If the data is similar, then column binding is more straingt forward, fast and less memory intensive. In this case, create one unique column that has all the enteries that will be used for merging, 
# and merge that column to all the dataframes. This can create a consistent column-entry for easy cbind() implementation.


# Another way to do multi-merge;
# Source - https://stackoverflow.com/a/8096127
reshape::merge_all(list_of_dataframes, ...)
reshape::merge_recurse(my.list)
# Speed comparison: the following test shows that merge_all() is the most efficient technique.

#==================================
# Negation, it should be the part of R-base
#=================================
# to negate the function in r 
"%notin%" <- function(x,table) match(x,table, nomatch = 0) == 0
'%notin%' <- Negate('%in%')

#==================================
# Change warning behavior
#=================================

# change the bahavior of warnings
# It may be useful to specify options(warn=2, error=recover)
# As mentioned by @plannapus, warn=2 will upgrade warnings to errors; error=recover will drop you into a debug/browser mode at the point where the warning (now upgraded to an error) occurred. 
# (Use options(warn=0, error=NULL) to restore the original settings.)

#==================================
# Trim white spaces 
#=================================
# A simple function to remove leading and trailing whitespace:
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

x <- c(" lead", "trail ", " both ", " both and middle ", " _special")
## gsub function from https://stackoverflow.com/a/2261149/7941188 
## this is NOT the function from user Bernhard Kausler, which uses 
## a much less concise regex 
gsub_trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#==================================
# String formating and editing
#=================================
# padding a string
des <- stringr::str_pad(des, 4, pad = "0")

#replace multiple patterns in name column
df$name <- gsub('A', 'Andy',
           gsub('B', 'Bob',
           gsub('C', 'Chad', df$name)))

ifelse(match(des_geoa$des, colnames(mat)), 
       gsub(., des_geoa$geoa, .), colnames(mat) )
des2 <- des[match(des, colnames(mat))]

#==================================
# GEO related, and other expression sets related
#=================================
gse <- getGEO("GSE33814", GSEMatrix = TRUE, 
                destdir="E:/path_to_dir/geo_data",
                getGPL = FALSE)

if (length(gse) > 1) 
	idx <- grep("GPL570", attr(gse, "names")) 
						else idx <- 1
gse <- gse[[idx]]

# custom CDF install; may be relevant "https://www.biostars.org/p/67294/"
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgcdf_22.0.0.tar.gz", type="source", repos=NULL)
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgprobe_22.0.0.tar.gz", type="source", repos=NULL)
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezg.db_22.0.0.tar.gz", type="source", repos=NULL)

# convert series matrix (expression set) to annotatedDataFrame or ExpressionFeatureSet
tADF <- AnnotatedDataFrame(#data = exprs(gse109),
                   pData(phenoData(gse109)))
                   
#new("ExpressionFeatureSet", exprs= exprs(gse109))
gse109_efs <- new("ExpressionFeatureSet",
                  assayData = assayData(gse109),
                  exprs= exprs(gse109), 
                  phenoData = phenoData(gse109))#, 
                  #featureData = featureData(gse109), 
                  #experimentData = experimentData(gse109), 
                  #annotation = gse109@annotation)#,
                  #platform = gse109@annotation)

#stopifnot(validObject(CelData))
validObject(gse109)


# create elist object to analyse Illumina microarray expression data; while analysing GSE109211, however, not all the variables were available at that time.
gse190_elist <- new("EListRaw")

gse190_elist@.Data[[1]] <- 'illumina'
?gse190_elist@.Data[[2]] <- targetinfo
gse190_elist@.Data[[3]] <- wdgs[[1]]
gse190_elist@.Data[[4]] <- gse109_data_exp
gse190_elist@.Data[[5]] <- NULL
gse190_elist$E <- gse109_data_exp

gse190_elist$targets <- targetinfo
gse190_elist$genes <- wdgs[[1]]
gse190_elist$other$Detection <- gse109_data_pvalue


# reconstruct the CelData; FeatureExpressionSet
# The FeatureSet class is VIRTUAL. Therefore users are not able to create instances of such class. 
# Objects for FeatureSet-like classes can be created by calls of the form: 
# new(CLASSNAME, assayData, manufacturer, platform, exprs, phenoData, featureData, experimentData, annotation, ...). 
# But the preferred way is using parsers like read.celfiles and read.xysfiles.

keep_row_pdata <- grepl("_2", rownames(pData(CelData)), ignore.case = T)

pData(CelData) <- pData(CelData)[!keep_row_pdata, ]
rownames(pData(CelData))  <- sub(pattern = "_1",replacement = "", x = rownames(pData(CelData)) )

# expression feature set
efs <- new("ExpressionFeatureSet",
           manufacturer = "Affymetrix",
           exprs = as.matrix(celdata_exp),
           phenoData =   phenoData(CelData),
           featureData = featureData(CelData),
           annotation = "pd.hg.u133a")

eCelData <- efs
validObject(eCelData)

#==================================
# Affymatrix can be manipulated in multiple ways
#=================================
affyids <- gse37_deg$ID
library(hgu133plus2.db)
columns(hgu133plus2.db)

select(hgu133plus2.db, keys = gse37_deg$ID, columns = c("SYMBOL"), keytype = "PROBEID")

x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
test <- clusterProfiler::bitr(gse37$ID, fromType = "PROBEID", toType = "SYMBOL",
                      OrgDb = "hgu133plus2.db", drop =TRUE)

# for easier access, you can search ID type in the biomaRt object, 
x = biomaRt::listAttributes(ensembl)
x[grep("affy", x$description, ignore.case = T),]
						
mart <- biomaRt::useEnsembl("ensembl","hsapiens_gene_ensembl")
converted_ID <- biomaRt::getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol'), # data types you want.
      filters = 'affy_hg_u133_plus_2',  # external_gene_name, or type of IDs you have.
      values = gse$ID,           		# list of ID
      mart = mart)

#==================================
# Convert geneIDs from ENSEM to ENTREIDs
#==================================

# While converting the names, BiTr usually returns a df with two columns, "fromType" & "toType", and it is oftenly difficult to match them
# to the original dataset. This is a simple function around that concept that BiTr should bind to converted IDs back to the original df 
# for easier downstream analysis

bitr2 <- function(df) {
    stopifnot(class(df)  == "data.frame")
    message("Input is not dataframe")
    items <- row.names(df)
    itemsID <- clusterProfiler::bitr(items, fromType="ENSEMBL", 
                                     toType=c("ENTREZID", "SYMBOL"),
                                     OrgDb=organism, drop=TRUE)
    df <- merge(df, itemsID, by.x = 0, by.y="ENSEMBL", all.x = TRUE)
    df <- transform(df, log2FoldChange = as.numeric(log2FoldChange), 
                    ENTREZID = as.numeric(ENTREZID))
    df <- df[complete.cases(df),]
    df = df[!duplicated(df$ENTREZID),]
    message("after removing NA: ", dim(df))
    df
}

# another way converting probeIDs to other IDs
require(hgu133a.db)

annotMaster1 <- select(hgu133a.db,
  keys = keys(hgu133a.db, 'PROBEID'),
  column = c('PROBEID',  'SYMBOL',  'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')

#==================================
# Feature selction and machine learing
#==================================
# First, identify the highly correlated attributes to save time, generally, > abs(0.75) or higher
# attributes should be removed; the correlation can be computed using base function "cor"

corrMatrix <- cor(df) # output is all * all cor matrix of same dim

highlyCorrelated <- caret::findCorrelation(corrMatrix, cutoff = 0.75)

# remove those features that have > 0.75 correlation 
dresist_lung_corr75 <- dresist_lung[, -highlyCorrelated]

# createDataPartition() function from the caret package to split the original dataset into
# a training and testing set and split data into training (80%) and testing set (20%)

parts = createDataPartition(dresist_lung$sampleType, p = 0.8, list = F)
caret.train = dresist_lung_corr75[parts, ]
caret.test = dresist_lung_corr75[-parts, ]
x_train = caret.train[, -length(dresist_lung_corr75)]
y_train = caret.train[,  length(dresist_lung_corr75)]

# specifying the CV technique as well as the random forest algorithm which will be passed into 
# the recursive feature eleminiation (rfe) rfe() function in feature selection

control_rfe = rfeControl(functions = rfFuncs, 	# random forest
                         method = "repeatedcv", # repeated cv
                         repeats = 10, 		# number of repeats
                         number =  10) 		# number of folds

# Performing RFE
result_rfe = rfe(x = x_train, 
                 y = as.factor(y_train), 
                 sizes = c(1:length(x_train)),
                 rfeControl = control_rfe)

# summarising the results
result_rfe

# creating a model using these features
fitControl <- trainControl(## 10-fold CV
    method = "repeatedcv",
    number = 10,
    repeats = 10,
    search = "random")  # hyper-parameters random search 

model.cv <- train(sampleType ~ cg02855924+cg26325335+cg27071460+cg14823851+cg01573747,
                  data = caret.train,
                  method = "glmnet", family = "binomial",
                  trControl = fitControl,
                  preProcess = c('scale', 'center'),
                  na.action = na.omit)

# print the model
print(model.cv)

plot(model.cv)

# score for TEST data set
class <- predict(model.cv, caret.test)
probs <- predict(model.cv, caret.test, 'prob')

# bind wtih actual data for easier inspection
TEST.scored <- cbind(caret.test,class,probs) %>% mutate(data = "TEST")

#==========================
# Linear modeling related
#==========================
# significant difference between models; likelyhood ratio test
anova(model.cv.f5.glm, model.cv.f4.glm, test = "Chisq")
anova(model.cv.f3.glm, model.cv.f5.glm, test = "Chisq")
anova(model.cv.f3.glm, model.cv.f4.glm, test = "Chisq")

lmtest::lrtest(model.cv.f3.glm, model.cv.f1.glm)

# pseudo R-squred R2 test; Most notable is McFadden’s R2, which is defined as 1−[ln(LM)/ln(L0)]
# where ln(LM) is the log likelihood value for the fitted model and ln(L0) is the 
# log likelihood for the null model with only an intercept as a predictor. 
# The measure ranges from 0 to just under 1, with values closer to zero indicating
# that the model has no predictive power.

pscl::pR2(model.cv.f5) # look for McFadden; if the model trained with glm and not with caret

# Statistical Tests for Individual Predictors
# Wald Test
# A wald test is used to evaluate the statistical significance of each coefficient
# in the model and is calculated by taking the ratio of the square of the
# regression coefficient to the square of the standard error of the coefficient. 
# The idea is to test the hypothesis that the coefficient of an independent variable in 
# the model is significantly different from zero. If the test fails to reject the null hypothesis, 
# this suggests that removing the variable from the model will not substantially 
# harm the fit of that model.

survey::regTermTest(model.cv.f5.glm, "ENSG00000169398")
survey::regTermTest(model.cv.f5.glm, feat5[5])
survey::regTermTest(model.cv.f3.glm, "ENSG00000169398")
survey::regTermTest(model.cv.f1.glm, "ENSG00000169398")



pred <- predict(model.cv.f5, newdata = caret.test)
accuracy <- table(pred, caret.test[, "status2"])
sum(diag(accuracy))/sum(accuracy)

confusionMatrix(data=pred, caret.test$status2)# not working


# ROC for single variable
f1 <- pROC::roc(status2 ~ model.cv.f5, data = caret.test) 
plot(f1)

#==================================
# Model evaluation methods
#==================================

# for the entire model

auc_check <- function(testModel, 
                      testData,
                      type = "prob",
                      targetColumn = "status2",
                      plot = TRUE) {
    require(ROCR)
    
    prob = predict(testModel, newdata = testData, type = type)[,2]
    
    pred = ROCR::prediction(prob, as.data.frame(testData[ , targetColumn]))
    perf = ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    if (plot) {
        plot(perf)}
    auc = ROCR::performance(pred, measure = "auc")
    auc = auc@y.values[[1]]
    return(auc) 
    
}

auc_check(testModel = model.cv.f5, testData = caret.train, plot = FALSE)


auc_check2 <- function(model, 
                       testData, 
                       trainData, 
                       targetClass = "targetCol")
{
    
    require(dplyr)
    require(yardstick)
    class <- predict(model, testData)
    probs <- predict(model, testData, "prob")
    
    TEST.scored <- cbind(testData, class, probs) %>% 
        mutate(data = "TEST")
    
    # score TRAIN
    class = predict(model, trainData)
    probs = predict(model, trainData, 'prob')
    
    TRAIN.scored = cbind(trainData, class, probs) %>% 
        mutate(data = "TRAIN")
    
    TRAIN_TEST_scored = rbind(TRAIN.scored, TEST.scored)
    TRAIN_TEST_scored[, targetClass] = as.factor(TRAIN_TEST_scored[, targetClass])
    
    #library(fpps)
    TRAIN_TEST_scored %>%
        group_by(data) %>%
        roc_auc(truth = targetClass, "Rec")
}

auc_check2(model.cv.f5, caret.test, caret.train, targetClass = "status2")


# confusion matrix will not work if type != "raw", 

pred.prob <- predict(model.cv.f5, newdata = caret.test, type = "raw")
caret::confusionMatrix(data = pred.prob, caret.test$status2, positive = "Rec")

#==================================
# ROC curve in R
#==================================
library(pROC)

roc.cg00074348 <- pROC::roc(data = betas.fs.2, 
		response = "sampleType",
		predictor = "cg00074348",
		ret = c("roc", "coords", "all_coords"),
		ci = TRUE, plot = TRUE)

pROC::plot.roc(roc.cg00049664,
	xlim=if(roc.cg00049664$percent){c(100, 0)} else{c(1, 0)},
	ylim=if(roc.cg00049664$percent){c(0, 100)} else{c(0, 1)})

roc.cg <- pROC::roc(data = test, #betas2.fs, 
		response = "sampleType",
		predictor = "y_pred",
		ret = c("roc", "coords", "all_coords"),
		ci = TRUE, plot = TRUE)

pROC::plot.roc(roc.cg,
	xlim=if(roc.cg$percent){c(100, 0)} else{c(0, 1)},
	ylim=if(roc.cg$percent){c(100, 0)} else{c(0, 1)})


roc_df <- data.frame(
  TPR=rev(roc.cg$sensitivities), 
  FPR=rev(1 - roc.cg$specificities), 
  labels=roc.cg$response, 
  scores=roc.cg$predictor)						

#====================================
# Brier Score calculation
#====================================

calcBrierScore <- function(model, 
                           testdata = "testData",
                           targetCol = "status2")
    {
    
    require(DescTools)
    # first calculate the probabilities of prediction
    pred.prob = stats::predict(model, newdata = testdata, type = "prob")
    
    bScore1 = DescTools::BrierScore(as.numeric(testdata[ , targetCol]) -1, pred.prob[, 2])
    cat("From DescTools::BrierScore: ", bScore1)
    
    # Brier Score is similar to MSE; therefore
    f_t = pred.prob[,2]
    o_t = as.numeric(testdata[ , targetCol ] ) - 1 # species are the target class
    bScore2 = mean((f_t - o_t)^2)
    cat("\nFrom formula: ", bScore2)
    }

calcBrierScore(model.cv.f5, testdata = caret.test)

# https://stackoverflow.com/questions/61014688/r-caret-package-brier-score
# I use the Brier score to tune my models in caret for binary classification. 
# I ensure that the "positive" class is the second class, which is the default when 
# you label your response "0:1". Then I created this master summary function, 
# based on caret's own suite of summary functions, to return all the metrics I want to see:

BigSummary <- function (data, lev = NULL, model = NULL) {
    pr_auc <- try(MLmetrics::PRAUC(data[, lev[2]],
                                   ifelse(data$obs == lev[2], 1, 0)),
                  silent = TRUE)
    brscore <- try(mean((data[, lev[2]] - ifelse(data$obs == lev[2], 1, 0)) ^ 2),
                   silent = TRUE)
    rocObject <- try(pROC::roc(ifelse(data$obs == lev[2], 1, 0), data[, lev[2]],
                               direction = "<", quiet = TRUE), silent = TRUE)
    if (inherits(pr_auc, "try-error")) pr_auc <- NA
    if (inherits(brscore, "try-error")) brscore <- NA
    
    rocAUC <- if (inherits(rocObject, "try-error")) {
        NA
    } else {
        rocObject$auc
    }
    tmp <- unlist(e1071::classAgreement(table(data$obs,
                                              data$pred)))[c("diag", "kappa")]
    out <- c(Acc = tmp[[1]],
             Kappa = tmp[[2]],
             AUCROC = rocAUC,
             AUCPR = pr_auc,
             Brier = brscore,
             Precision = caret:::precision.default(data = data$pred,
                                                   reference = data$obs,
                                                   relevant = lev[2]),
             Recall = caret:::recall.default(data = data$pred,
                                             reference = data$obs,
                                             relevant = lev[2]),
             F = caret:::F_meas.default(data = data$pred, reference = data$obs,
                                        relevant = lev[2]))
    out
}

# Now I can simply pass  summaryFunction = BigSummary in trainControl and then 
# metric = "Brier", maximize = FALSE in the train call.
						
#==================================
# boxplot
#==================================
checkDiffExp <- function(geneName, df) {
    normal = df[rownames(df) == geneName, samples.solid.tissue.normal]
    tumor  = df[rownames(df) == geneName, samples.primary.tumour]
    y_max = ifelse(any( max(normal) | max(tumor)) > 25000, 25000, max(max(normal),  max(tumor))) 
    boxplot(normal, tumor, ylim = c(0, y_max),
            main= geneName,
            names = c("Normal", "Tumor"),
            xlab= "Type",
            ylab= "Expression Level") }

checkDiffExp("ENSG00000000938", dff)

#==================================
# cbind with diff number of rows
#==================================
# cbind requires same number of rows, if not it will fail. This is to avoid that, and it will fill in with NA
cbind.fill <- function(...) {
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
							}

#==================================
# package version and unloading of package. 
#==================================						  

# Unload name space, remove package from environment, unload
unloadNamespace("TCGAbiolinks")
detach("package:TCGAbiolinks", unload = TRUE, force = TRUE)

packageVersion("TCGAbiolinks")

#==================================
# Converting UUID values to file names, 
#==================================

manifest <- read.csv("TCGA_data/Methylation_betaValues/gdc_manifest.2024-03-19.txt",
                            header = T, sep = "\t")

manifest_length= nrow(transc.manifest)
id= toString(sprintf('"%s"', transc.manifest$id))

Part1= '{"filters":{"op":"in","content":{"field":"files.file_id","value":[ '
Part2= '] }},"format":"TSV","fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id","size":'
Part3= paste(shQuote(manifest_length),"}",sep="")
Sentence= paste(Part1,id,Part2,Part3, collapse=" ")
write.table(Sentence,"Payload_betaVal.txt",quote=F,col.names=F,row.names=F)

curl --request POST --header "Content-Type: application/json" --data @Payload_betaVal.txt "https://api.gdc.cancer.gov/files" > File_metadata.txt

# check the intersection or differnces between two data set
intersect(rownames(col.data.man),  colnames(g12d_rna_seq_man2))
setdiff(rownames(col.data.man),  colnames(g12d_rna_seq_man2))


# to remove list with lower number of cpgs, 
lower_bval <- list()

for (i in 1:108) {
    if (dim(g12d_meth[[i]])[[1]] == 27578) {
        lower_bval <- append(lower_bval, g12d_meth[[i]][1, 1:2])
    }
}

#==================================
# migrate R libraries, move R libraries, install R libraries
#==================================
# In original installation, get the non-default package list:
save.pkg.list <- installed.packages()[is.na(installed.packages()[ , "Priority"]), 1]
save(save.pkg.list, file="pkglist.Rdata")
# If you want to use remove.packages() at this point it's fine. 
# Or just delete their directories.

load("pkglist.Rdata")
install.packages(save.pkg.list)

#==================================
# InfiniumMethylation lib to convert probe id to gene 
#==================================

library(FDb.InfiniumMethylation.hg19)
# list all the contents in the package
ls('package:FDb.InfiniumMethylation.hg19')
 [1] "FDb.InfiniumMethylation.hg19" "get27k"                      
 [3] "get450k"                      "getNearest"                  
 [5] "getNearestGene"               "getNearestTranscript"        
 [7] "getNearestTSS"                "getPlatform"                 
 [9] "hm27ToHg19"                   "lift27kToHg19"     

hm450 <- get450k()
hm450k.probe2gene <- hm450[res$Probe_ID]
res.tss <- getNearestTranscript(hm450k.probe2gene)
res2.annotated.sig.tss2 <- getNearestTSS(res2.annotated.sig.gene)

#==================================
# GEDI data integration, and batch correction
#==================================						  
# Ref: https://github.com/MNStokholm/GEDI
library(GEDI)

dataFolders <- c("lung_data/",
                 "datasets/GSE_RAW/",
                 "datasets/GSE_RAW2/")

sources <- c("RNAseq", "affy", "affy")

# Read the data
PATH_TO_DATA_FOLDERS <- "E:/drug_resistance/"
datasets <- ReadGE(dataFolders, sources, path = PATH_TO_DATA_FOLDERS)
attr <- c("ensembl_gene_id", "affy_hg_u133a_2" , "affy_hg_focus")

# The datasets are integrated. The species Bos taurus is used
# Due to bug in dbplyr, BiomaRT is not working at this moment.

dat <- GEDI(datasets, attributes = attr, BioMart = TRUE,
            species = "hsapiens", path = PATH_TO_DATA_FOLDERS)

#==================================
# Linear/logistic regression analysis For feature selection
#==================================
#probes should be in column, samples in rows, and additional column for sample type
library(FSinR) # feature selection in R

# for hybrid search method, two evlaution methods are used:
# (1) filter method
# (2) wrapper method
#
# at the time of writing, only one algorithm was supported to perform hybrid search
# hybrid_search_method <- hybridSearchAlgorithm('LCC')

# Generates the first filter evaluation function (individual or set measure)
f_evaluator <- filterEvaluator('determinationCoefficient')

# Generates the second wrapper evaluation function (mandatory set measure)
resamplingParams <- list(method = "cv", number = 10)
fittingParams 	 <- list(preProc = c("center", "scale"), metric="Accuracy", tuneGrid = expand.grid(k = c(1:20)))
w_evaluator <- wrapperEvaluator("knn", resamplingParams, fittingParams)

# Generates the hybrid search function
hybrid_search_method <- hybridSearchAlgorithm('LCC')

# Runs the hybrid feature selection process
featureSel <- hybridFeatureSelection(df, "sampleType", hybrid_search_method, f_evaluator, w_evaluator)

roc.curve <- pROC::roc(	data = df, 
		response = "sampleType",
		predictor = "cg12483545",
		ret = c("roc", "coords", "all_coords"),
		ci = TRUE, plot = TRUE)

plot.roc(roc.curve,
	xlim=if(roc.curve$percent){c(100, 0)} else{c(1, 0)},
	ylim=if(roc.curve$percent){c(0, 100)} else{c(0, 1)})

roc_df <- data.frame(
  TPR=rev(roc.cg$sensitivities), 
  FPR=rev(1 - roc.cg$specificities), 
  labels=roc.cg$response, 
  scores=roc.cg$predictor)

#==================================
# Plot various aspects in linear/log regression
#==================================

# needs to updat
						  
#==================================
# Machine learning in R
#==================================
library(mlbench)
library(caret)

#================================== 
# detach a package in R
#==================================
detach("package:pscl", unload=TRUE)
#You can also use the unloadNamespace command, as in:
unloadNamespace("pscl")

#================================== 
# ggplot in R
#==================================
library(reshape2)
library(dplyr)

# incase label overlapps in ggplot and ggrepel
options(ggrepel.max.overlaps = 10)

dds_normCount2[, c(2:21)] %>% 
    filter(row.names(dds_normCount2) == "ENSG00000000170") %>%
    melt() %>% mutate("type" = c(rep("normal", 10), rep("obse", 10))) %>%
    ggplot(aes(x = type, y = log2(value), fill = type)) + 
    geom_errorbar(stat = "summary", width = 0.1, color = "black", alpha = 1.5) +
    stat_summary(geom = "point", fun = mean, color = "black") +
    #geom_point(position = position_jitter(width = 0.1), shape = 18, size = 4) +
    scale_color_brewer(palette = "Set2") +
    #geom_boxplot() + 
    geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=1.2) + 
    theme_classic() + 
    labs(title="PPARg",
         x ="Type", 
		 y = "Normalized expression")

# save figure; can auto save the last generated plot
ggsave(filename = "pparg_from_gse162653.pdf", 
	   plot = get_last_plot()) # or the plot name

# plot a barplot for top 10 enriched terms ordered by q-values
ggplot(cpGO_df[1:10, ], aes(x = -log10(qvalue[1:10]), y =  reorder(Description[1:10], -log10(qvalue[1:10]) )) + 
	geom_bar(stat = "identity") + 
	theme_classic()

# create a star bust plot 
ggplot(data = deg_dm, aes(x = Est_groupTumor, y = log2FoldChange)) + 
    geom_point(colour = "white", shape = 21, size = 2, aes(fill = factor(dataType))) + 
    geom_hline(yintercept = c(1, -1), linetype = "dashed", color = "red") + 
    geom_vline(xintercept = c(0.3, -0.3), linetype = "dashed", color = "blue") +
    scale_x_continuous(breaks = c(0, 0.3, -0.3), minor_breaks=NULL) +
    scale_y_continuous(breaks = c(seq(0,8,2), seq(0,-8,-2)), minor_breaks=NULL) +
    xlab("Methylation differences") + 
	   ylab("log2FoldChange") +
    theme_bw()  	   

# ggplot; volcano plot in ggplot2
ggplot(data = deg_dm, aes(x = Est_groupTumor, y = log2FoldChange)) + 
    geom_point(colour = "white", shape = 21, size = 2, aes(fill = factor(dataType))) + 
    theme_classic() + 
    geom_hline(yintercept = c(1, -1), linetype = "dashed", color = "red") + 
    geom_vline(xintercept = c(0.3, -0.3), linetype = "dashed", color = "blue") +
    scale_x_continuous(breaks = c(0, 0.3, -0.3), minor_breaks=NULL) +
    scale_y_continuous(breaks = c(seq(0,8,2), seq(0,-8,-2)), minor_breaks=NULL) +
    xlab("Methylation differences") + ylab("log2FoldChange")
	   
	geom_boxplot() +  				# creating box plot
    geom_dotplot(binaxis='y', stackdir='center', 	# craeting dot plot
                     stackratio=1.5, dotsize=1.2) + 

# create plots from files
plist <- list()
for (i in list.files("results_and_fingures/", 
                     pattern = "GO_Biological", 
                     recursive = T, 
                     full.names = T)) {
    
    title = stringr::str_split_i(i, pattern = "/", 2)
    
    if( title == "Common_genes_between_Reactome_KEGG" ) {
        stringr::str_split_i(title, pattern = "_", 1)
    } else { }
    
    print(title) # it should be cat()
    first <- read.csv(i, sep = "\t", header = T, stringsAsFactors = F)[, c(1,4)]

    plist[[i]] <- ggplot(first[1:10, ], aes(x = -log10(Adjusted.P.value[1:10]), y =  reorder(Term[1:10], -log10(Adjusted.P.value[1:10])))) + 
        geom_bar(stat = "identity", fill = "coral3") + 
        theme_classic() 			+ 
        xlab("-log10(adj.p-value")	+ 
        ylab("Term") 				+ 
        scale_y_discrete(labels = function(x) {
            is_long <- nchar(x) > 25
            x[is_long] <- paste0(substr(x[is_long], 1, 25), ".") # this can trim the label, from very long to a reasonable length
            x }) +
        ggtitle(title)
}

#===================================================
# EnrichR for RNA seq data
#===================================================

#install.packages("enrichR")
library(enrichR)

# by default, species == human 
# Then find the list of all available databases from Enrichr.
# dbs <- listEnrichrDbs()
# listEnrichrSites()

enrichr_GO <- function(geneNames, fName) {
    
    websiteLive <- TRUE
    
    dbs <- c("GO_Biological_Process_2021",	# 3 
             "KEGG_2021_Human",
             "MSigDB_Hallmark_2020"			# 4
    )
    
    if (websiteLive) { enriched <- enrichr(geneNames, dbs) }
    
    for (i in 1:length(dbs)){
        
        write.table(enriched[[i]][which(enriched[[i]]$Adjusted.P.value < 0.05),], paste0(fName,"_",dbs[i], ".txt", sep=""), sep = "\t")
    }
     
    message("\n Saving files in current working directory: ", getwd(), "\n")
    # For plot the GO enrichments
    #if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
    #enriched
    }

plot_from_file <- function(pathDir, fname){
    require(ggplot2)
    
    file = read.table(paste0(pathDir, fname), sep = "\t", header = T, 
                      stringsAsFactors = F)
    file = file[ order(file[, c(4)], decreasing = F), ]
    print(head(file))
    ggplot(file[1:10, ], aes(x = -log(Adjusted.P.value), y = reorder(Term, -log10(Adjusted.P.value)))) +
        geom_bar(stat = "identity", fill = "orange", color = "orange") +
        theme_classic()  
        
    ggsave(filename = paste0(pathDir, fname, ".pdf"))
}

pathDir <- "result/rna_seq_enrichment_EnrichR/"
plot_from_file(pathDir, "normalUp_MSigDB_Hallmark_2020.txt")
plot_from_file(pathDir, "cancerUp_MSigDB_Hallmark_2020.txt")

#================================== 
# Heatmap in R
#==================================

# for heatmap; need to create a comprehensive structure
heatmap ==> default or from base

pheatmap::pheatmap(exprs(Normalize2), 
			annotation_col = ppp_p[, "response", drop = F] # This is the column/row annotation.
			)

ComplexHeatmap::Heatmap(as.matrix(hm.all2[1:10,]),
                        na_col = "grey",
                        column_dend_side = "top",
                        col = circlize::colorRamp2(c(-3, 0, 3), c("Darkblue", "white", "red")),
                        heatmap_legend_param = list(color_bar = "continuous"))

#==================================
# nested ifelse
#==================================

# ifelse() can be nested in many ways:

ifelse(<condition>, <yes>, 
	   ifelse(<condition>, <yes>, <no>)
		)

ifelse(<condition>, ifelse(<condition>, <yes>, <no>), <no>)

ifelse(<condition>, 
       ifelse(<condition>, <yes>, <no>), 
       ifelse(<condition>, <yes>, <no>)
      )

ifelse(<condition>, <yes>, 
       ifelse(<condition>, <yes>, 
              ifelse(<condition>, <yes>, <no>)
             ) )

#==================================
# for Geometric Mean calculation
#==================================
exp(mean(log(x))) 

# this and the foloowing function is same or in function from
geometric.mean <- function(x, na.rm=TRUE) { 
		exp(mean(log(x),na.rm=na.rm))		}

#==================================
# Install new fonts in R
#==================================
# to install new fonts in R: have to done all three steps
library(extrafont)
font_import()
loadfonts(device = "win")


#==================================
# PCA
#==================================

# PCA and plot
# quality check for normalized data
  mdaPcaRma <- Biobase::exprs(mdaDataRma)
  
  mdaPcaRma <- prcomp(t(mdaPcaRma), scale. = T)
  
  PercentVarRma <- round(100 * mdaPcaRma$sdev^2 / sum(mdaPcaRma$sdev^2), 1)
  SdRatioRma <- sqrt(PercentVarRma[2] / PercentVarRma[1])
  
  mdaPcaRmaGG <- data.frame(PC1 = mdaPcaRma$x[, 1], 
                            PC2 = mdaPcaRma$x[, 2],
                            Disease = pData(mdaData)$`status of disease:ch1`,
                            Phenotype = pData(mdaData)$`response to therapy:ch1`,
                            Individual = pData(mdaData)$`geo_accession`,
                            Gender = pData(mdaData)$`gender:ch1`)
    
  
  ggplot(mdaPcaRmaGG, aes(PC1, PC2)) +
    geom_point(aes(shape = Disease, colour = Phenotype), size = 5) +
    ggtitle("PCA plot of the calibrated, summarized data") +
    xlab(paste0("PC1, VarExp: ", PercentVarRma[1], "%")) +
    ylab(paste0("PC2, VarExp: ", PercentVarRma[2], "%")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = SdRatioRma) +
    scale_shape_manual(values = c(1)) + 
    scale_color_manual(values = c("darkorange2", "dodgerblue4", "green")) +
    theme_classic()
  
#=========================================
# Survival analysis
#=========================================
library(survival)
coxph(Surv(time = days, 
           event = as.numeric(as.character(vital_status_0alive_1dead)), 
             type = "right") ~ gene1 + gene2 + gene1:gene2, data = survival_data) # gene1:gene2 is an interaction term
				   
# survival analysis parallel
library(survival)
library(RegParallel)

res <- RegParallel(
  data = vsd_survival,
  formula = 'Surv(vital_status_0alive_1dead, days_to_death_2) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(vsd_survival)[1:32056], #c("vital_status_0alive_1dead", "days_to_death_2"), #colnames(vsd_survival)[ncol(vsd_survival)-1:ncol(vsd_survival)],
  blocksize = 2000,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95)

res <- res[!is.na(res$P),]
res

#========================================
# helper function to select one column, while analyzing TCGA data
#========================================

# Function to select one column
keep_one_column <- function(input_df, term){
  
  col_list <- grep(term, colnames(colData(input_df)), ignore.case = T, value = T)
  print(col_list)

  input_df <- as.data.frame(colData(input_df)[, col_list])
  print(head(input_df))

  if (any(grepl("days", col_list))) {
                  mat_col <- input_df %>% # as.data.frame() %>%
                    mutate_if(is.character, as.numeric)
  } else {
    mat_col <- input_df
  }
  
    # if (length(col_list) > 1){
    # mat_col <- colData(input_df)[ col_list ]
    
    cat("\n\ncolumn names are:\n\n ", colnames(mat_col))
    mat_col$merged <- dplyr::coalesce(mat_col[, 1], mat_col[, 2] )
    
  return (mat_col)
}

#=========================================
# create a random dataset; for testing and evlaution purposes
#=========================================
create_random_data <- function(row_num, col_num, output = "matrix", range = c(0, 10) ) {
	# a small function to create test data set
	# the out put should be either matrix of df, so if you want it to be df, just write anything.
	# though, it is a crude way to do, but it works.
	# It should be used for matrix of df type of random data set, as for simple numeric or character vectors, 
	# base::sample can be used i.e 
	# base::sample(letters, 5) --> "q" "p" "h" "v" "u"
	# base::sample(0:10, 5) --> 3 5 6 1 0

# added a check
  stopifnot("rows can not be 0"  = row_num > 0 ,
            "cols can not be 0"  = col_num > 0)
	
 out <- replicate(col_num, sample(range[1]:range[2], size= row_num), simplify=FALSE)
 out <- do.call("cbind", out)
  
  if (output != "matrix"){
    out <- as.data.frame(out)
    # assign col names, in case of df, the default value is "V" with increasing number.
    column_names <- paste0(rep("col_", col_num), seq_len(col_num))
    names(out) <- column_names
  }
  
  return(out)
}

create_random_data(3, 5, output = "matrix") # it will create matrix
create_random_data(3, 5, output = "mat") # it will create data frame

# need to adjust
DF <- data.frame(A = rnorm(5), # normal distribution
                 B = rnorm(5),
                 C = rnorm(5),
                 D = rnorm(5),
                 E = rnorm(5),
                 F = rnorm(5))

# to create random distribution with lower bound (-0.68) and upper bound (0.82)
random_dist <- truncnorm::rtruncnorm(n=108, a=-0.68, b=0.82)#, mean=0.003, sd=0.292)				   

#=========================================
# Set col to row and delete it. 
#=========================================
setColToRow <- function(df = df, 
                        col.name = "col.name") {
    # if the column enteries are duplicated, it will make them uniq, otherwise, rownames can not be set.
    stopifnot("Column name is not present in data frame: " = col.name %in% colnames(df))
    
    # unique the names; if not unique, rowname can not be set
    df[[col.name]] = make.unique(df[[col.name]])
    
    if(!is.data.frame(df)) {
        stop("df is not cohercible to df")
    } else
    {
        rownames(df) = df[[col.name]]
        df[[col.name]] = NULL
        
        return(df)
    }
}




#=========================================
# add snippet to the RStudio
#=========================================
# To paste windows path to Rstudio, a snippet has been added in RStudio. 
# Under windos system, the snippet can be activated by shift-tab
# here, pp is snippet name, write pp and then press shift-tab will paste the altered path (if copied).
# copy (file path) --> pp --> shift-tab
# the snipped is as follow (it will substitute the forward slash with the bacslashes):
snippet pp
	"`r gsub('"', "", gsub("\\\\", "/", readClipboard()))`"

#(2)
snippet .
	%>% 
#(3)
snippet tv
	library(tidyverse)
# (4)
snippet ss
	#=========================================
	#
	#=========================================












