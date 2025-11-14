library(stringr)
library(dplyr)

# This is for general history and commands list

# select max value between two columns; pmax is parallel maximum
common1$log10pvalue <- pmax(common1$log10pvalue.x, common1$log10pvalue.y)

#Run this command in case many background files are open e.g. plots or pdf files
dev.off()

# to chek the type of individual coloumn types
str(data_frame)


samples <- colnames(lung.prot)
library(stringr)
samples <- gsub("(?<=QC1)[^;]*", "", samples, perl = TRUE)
samples2 <- samples[2:973]
samples2 <- grepl("Unshared", samples2)
sum(grepl("LateStageTumor_pool.QC1", samples))
sum(grepl("Unshared.", samples))
samples2 <- grep("Unshared.", samples, value =T)
samples2 <- str_remove(samples2, "Unshared.Log.QC2_LateStageTumor_pool.QC1")
samples2 <- samples2[nzchar(samples2)]


# to create random distribution with lower bound (-0.68) and upper bound (0.82)
random_dist <- truncnorm::rtruncnorm(n=108, a=-0.68, b=0.82)#, mean=0.003, sd=0.292)

# apply Wilcoxon test (this does not assume normal distribution)
wilcox.test(as.numeric(common2[1, 2:193]), random_dist, paired = F, alternative = "two.sided")

#===================================
# t-test rowwise (Welch t-test)
#===================================
# by default, R performs Welch Two Sample t-test
x$stat <- sapply(1:nrow(x), function(i) 
	t.test(	as.numeric(as.character(unlist(x[i,2:4]))), 
			as.numeric(as.character(unlist(x[i,8:10])))
			)[c("p.value")])

#===================================
# Limma for gene expression
#===================================
# limma for microarray data analysis
# limma is for continuous data while the DESeq2 and EdgeR are for count data
# Moreover, in microarray, the Robust Muti-array analysis (RMA) converts the intensities into log2 form, making it easy to find log2FC
# and simply running the t test, the significance can be calculated

design_limm <- model.matrix(~ factor(mapping3$sampleType))
fit <- lmFit(eset, design_limm)
ebayes <- eBayes(fit)

lod <- -log10(ebayes[["p.value"]][, 2])
mtstat <- ebayes[["t"]][, 2]

# top significant genes
tab <- topTable(ebayes, coef=2, adjust="fdr", n=10)

				 
# create expression set for GEOquery 
# it was a test, but the GEO data is very unpredictable/messsy
# making it really hard to convert it into ExpressionSet

eset <- ExpressionSet(assayData = as.matrix(mat.gse126848),
                      phenoData =  Biobase::AnnotatedDataFrame(phenoData(test2[[1]])))

phenoData <- new("AnnotatedDataFrame",
                 data=pData(test2[[1]])) 
                 #varMetadata=metadata)                      

exampleSet <- ExpressionSet(assayData=as.matrix(mat.gse126848),
                               phenoData=phenoData,
                               experimentData=experimentData(test2[[1]]),
                               annotation="hgu95av2")

# data.table fread is faster than regular read.csv
file <- data.table::fread("summed_tum_normal_refine_sam1.csv",
			sep = "\t", blank.lines.skip=TRUE, header = TRUE)

output1 <- Sys.glob("NCC_*_Proteome_KU_*\\OUTPUT\\")

# pattern based matching and retreiving the data
paths <- list.files(output1, 
					pattern= glob2rx("*summed_tum_normal_refine*.csv$*"),
                    full.names=T, recursive=T)

# Make a function to process each file, file name or identifier has to be appended to the respective columns. 
# after that, using sapply, this function can be used to all the files.
processFile <- function(f) {
  bname = strsplit(basename(f), '_')[[1]][1]
  df = data.table::fread(f, 
						 sep = "\t", 
						 blank.lines.skip=TRUE, 
						 header = TRUE, 
                         select = c(1,6:9), 
						 data.table = TRUE)
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

#===================================
# Merging
#===================================
# Multimerge, merge multiple dfs, if the number of dfs are large and the size is also big, it can run out of memory.
# allow.cartesian = TRUE is default in DF, while in data.table merge, it is FALSE, 
# you need to set it explicitly to proceed merging of data.tables.

# base R solution
merged_output <- Reduce(function(x, y) merge(x, y, all.x = TRUE, by = c("uniprot"), 
                                            allow.cartesian=TRUE), result[1:40])

# purr base solution
list_of_data %>% purrr::reduce(left_join, by = "row_names") # for purr based solution to merge multiple dataframes; but it needs a unique col name (e.g. row_names) in each df.

# data table approach,let see how efficient it is
library(data.table)
setkey(result_merge2, "uniprot")
df2_3 <- as.data.table(merge2)[as.data.table(merge1), on = "uniprot"]# allow.cartesian=TRUE ]

#result_merge <- result_merge[rowSums(is.na(result_merge[, 2:ncol(result_merge)])) == 0, ]
#result_merge <- result_merge %>%  dplyr::select(-starts_with("gene"))

lapply(result, dim)

#==================================
# Negation, it should be the part of R-base
#=================================
# to negate the function in r 
"%notin%" <- function(x,table) match(x,table, nomatch = 0) == 0
'%notin%' <- Negate('%in%')

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

# change the bahavior of warnings
# It may be useful to specify options(warn=2, error=recover)
# As mentioned by @plannapus, warn=2 will upgrade warnings to errors; error=recover will drop you into a debug/browser mode at the point where the warning (now upgraded to an error) occurred. 
# (Use options(warn=0, error=NULL) to restore the original settings.)


#==================================
# String formating and editing
#=================================
# padding a string
des <- stringr::str_pad(des, 4, pad = "0")

#replace multiple patterns in name column
df$name <- gsub('A', 'Andy',
           gsub('B', 'Bob',
           gsub('C', 'Chad', df$name)))

ifelse(match(des_geoa$des, colnames(mat.gse126848)), 
       gsub(., des_geoa$geoa, .), colnames(mat.gse126848) )
des2 <- des[match(des, colnames(mat.gse126848))]

#==================================
# GEO related
#=================================
gse <- getGEO("GSE33814", GSEMatrix = TRUE, 
                destdir="E:/NASH/geo_data_nash",
                getGPL = FALSE)

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
mart <-    biomaRt::useEnsembl("ensembl","hsapiens_gene_ensembl")
converted_ID <- biomaRt::getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol'), 
      filters = 'affy_hg_u133_plus_2',  #'external_gene_name', #
      values = gse37$ID,           # what you have
      mart = mart)


#==================================
# Feature selction and machine learing
#==================================
# First, identify the highly correlated attributes to save time, generally, > abs(0.75) or higher
# attributes should be removed; the correlation can be computed using base function "cor"

corrMatrix <- cor(df)

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

model.cv

plot(model.cv)

# score for TEST data set
class <- predict(model.cv, caret.test)
probs <- predict(model.cv, caret.test, 'prob')

# bind wtih actual data for easier inspection
TEST.scored <- cbind(caret.test,class,probs) %>% mutate(data = "TEST")

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
# biomaRt
#==================================
library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id"), 
                 values = probs,
                 mart = ensembl)


x = listAttributes(ensembl)
x[grep("affy", x$description, ignore.case = T),]

# mdata <- mdata[mdata$samp %in% colnames(cData2) , ]
smapT$status2 <- ifelse(grepl("Rec", smapT$status, fixed = T), "Rec", "Tum")

smapT$status2 <- ifelse(grepl("Non-Recurrent-NSCLC", smapT$status, fixed = T), "Tum", smapT$status2)    
smapT$status <- NULL


# unload name space, remove package from environment, unload
unloadNamespace("TCGAbiolinks")
detach("package:TCGAbiolinks", unload = TRUE, force = TRUE)

packageVersion("TCGAbiolinks")



reader.maf.files <- function(fpath, numb){
    bname = strsplit(basename(fpath), split = ".", fixed = TRUE,)[[1]][1]
    file = read.table(fpath, sep = "\t", header = TRUE)
    file = file[col.to.keep]
    colnames(file)[1] <- bname
    colnames(file)[3:length(file)] = paste0(colnames(file)[3:length(file)], "_", counter)
    counter <<- counter + 1
    return(file)
}

counter = 1
maf.luad <- lapply(files, reader.maf.files)


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


# ggplot ggplot2 
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
# 
#==================================
library(mlbench)
library(caret)

#================================== 
# list function in R package
#==================================
ls("package:enrichR")



intersect(Up, Dn)
# "Metazoa_SRP" "Y_RNA" 

#================================== 
# detach a package in R
#==================================
detach("package:pscl", unload=TRUE)
#You can also use the unloadNamespace command, as in:
unloadNamespace("pscl")


# reading excel files
readxl::read_excel

#================================== 
# ggplot in R
#==================================
library(reshape2)
library(dplyr)

dds_normCount2[, c(2:21)] %>% 
    filter(row.names(dds_normCount2) == "ENSG00000132170") %>%
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
# save figure
ggsave(filename = "pparg_from_gse162653.pdf", plot = )

#plot a barplot for top 10 enriched terms ordered by q-values
ggplot(cpGO_df[1:10, ], aes(x = -log10(qvalue[1:10]), y =  reorder(Description[1:10], -log10(qvalue[1:10]) )) + 
	geom_bar(stat = "identity") + 
	theme_classic()

# incase label overlapps in ggplot and ggrepel
options(ggrepel.max.overlaps = 10)

# In case, log2fc values of 2 diff gene are same, small increase small increment
#it will add 0.00001 to the second gene of the same log2fc value
mydata <- mydata %>%
		group_by(log2fc) %>% 
		mutate(log2fc2 = log2fc + seq(0, by=0.00001, length.out= n()))

#================================== 
# Heatmap in R
#==================================

# for heatmap
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
#
#==================================
# convert character to numeric in data frame
raw.ad[] <- sapply(raw.ad, as.numeric) ; it preserves both col names and row names, plus it requires [] on left side


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


# string split or you can remove the remaining part from ENSEMBL name, the version of ensembl id such as ENSG0000000000012.3. Remove .3, otherwise megeing will be difficult.
rownames(df) <- sub("\\..*", "", rownames(df))

strsplit use case
If you need to extract the first (or nth) entry from each split, use:
word <- c('apple-orange-strawberry','chocolate')
sapply(strsplit(word,"-"), `[`, 1)
#[1] "apple"     "chocolate"

Or faster and more explictly:

vapply(strsplit(word,"-"), `[`, 1, FUN.VALUE=character(1))
#[1] "apple"     "chocolate"

Both bits of code will cope well with selecting whichever value in the split list, and will deal with cases that are outside the range:

vapply(strsplit(word,"-"), `[`, 2, FUN.VALUE=character(1))
#[1] "orange" NA  


#drop columns that have 1 factor level
df[sapply(df, nlevels) > 1]

# for Geometric Mean calculation

exp(mean(log(x))) # this and the foloowing function is same

"geometric.mean" <- 
function(x, na.rm=TRUE) { 
		exp(mean(log(x),na.rm=na.rm))
						}


# compare all files in the list
outer(allFiles, allFiles, Vectorize(all.equal))

allFiles <- lapply(listOfFiles, function(x) readr::read_tsv(x,
                                                col_names = T,
                                                skip_empty_rows = T,
                                                trim_ws = TRUE)


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
  


# custom CDF install
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgcdf_22.0.0.tar.gz", type="source", repos=NULL)
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezgprobe_22.0.0.tar.gz", type="source", repos=NULL)
# install.packages("http://mbni.org/customcdf/22.0.0/entrezg.download/hgu133ahsentrezg.db_22.0.0.tar.gz", type="source", repos=NULL)






# get the expression set from the GEO dataset
gse109211 <- getGEO("GSE109211", GSEMatrix = T)

if (length(gse109211) > 1) idx <- grep("GPL570", attr(gse109211, "names")) else idx <- 1
gse109211 <- gse109211[[idx]]


# write.table(biomaRt::listAttributes(mart), file = "biomaRt_attributes_list.txt",
#             sep = "\t", quote = F)

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


# version of R
version

# check where lib are installed
.libPaths()



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


#=========================================
# Survival analysis
#=========================================
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
