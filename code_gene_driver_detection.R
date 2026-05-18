#######################################
#Codigo para driver detection en hepatocarcinoma
######################################
#Librerias
#install.packages("devtools")
library(devtools)
#install_github("im3sanger/dndscv")

#Carga de los datos
library(dndscv)
muts = read.table("~/LCG/2ndyear/3rd semester/Genomica funcional/data_mutations_maf.txt", header=T, sep="\t", stringsAsFactors=F)
head(muts)

#Dimensiones y numero de samples en el dataset
cat("The number of unique samples in the dataset is ", length(unique(muts$Tumor_Sample_Barcode)), ".") #Check number of unique samples
cat("\nThe number of mutations in the dataset is ", nrow(muts), ".") #Check number of mutations

#Grafico de samples con su numero de mutaciones
barplot(sort(table(muts$Tumor_Sample_Barcode)), ylab="Number of mutations",
        xlab="Donors", las=2, names.arg="")
abline(h = 500, col = "red", lty = 2)

#FIltra e identifica a los hipermutantes
number_muts_IDs <- table(muts$Tumor_Sample_Barcode) #We count the number of mutations of each donor
hypermuts <- muts[muts$Tumor_Sample_Barcode %in% names(number_muts_IDs[number_muts_IDs >= 500]), ] #We reduce the original dataset to those that are hypermutators
print(number_muts_IDs[unique(hypermuts$Tumor_Sample_Barcode)]) #We print the names and number of mutations of each hypermutator


#Dataframe del input de la prueba de abajo
input <- data.frame(
  sampleID = muts$Tumor_Sample_Barcode,
  chr  = muts$Chromosome,
  pos  = muts$Start_Position,
  ref  = muts$Reference_Allele,
  mut  = muts$Tumor_Seq_Allele2
)

#Prueba estadistica para detectar driver genes
dout = dndscv(input, max_muts_per_gene_per_sample=3,max_coding_muts_per_sample=500, outmats=T)

#Genes significativos
data_genes = dout$sel_cv[which(dout$sel_cv$qglobal_cv<0.1),]
head(data_genes)

#Genes bajo seleccion negativa
print(
  data_genes$gene_name[
    (data_genes$wmis_cv > 1 & data_genes$pmis_cv > 0.05) |
      (data_genes$wnon_cv > 1 & data_genes$pnon_cv > 0.05) |
      (data_genes$wspl_cv > 1 & data_genes$pspl_cv > 0.05) |
      (data_genes$wind_cv > 1 & data_genes$pind_cv > 0.05)
  ]
)

#Oncogenes
oncogenes = data_genes[data_genes$gene_name %in% c("CTNNB1", "NFE2L2"), ]
oncogenes

#Tumor supresors
tumor_s = data_genes[data_genes$gene_name %in% c("TP53", "AXIN1", "ARID1A", "KEAP1", "TSC2", "BAP1", "RB1", "ARID2", "ACVR2A"), ]

#Not defined
not_def = data_genes[data_genes$gene_name %in% c("ALB", "NLGN1", "IRX1", "CRIP3", "KRTAP3-3","KRTAP22-1"), ]

#Selected mutations 
selected <- function(coef, num){
  Prop = (coef-1)/coef
  cat("La proporcion de mutaciones en seleccion es de:" ,Prop, "\n")
  Mutsel = Prop*num
  cat("El numero de mutaciones seleccionadas es de: ", Mutsel, "\n")
}

##Para CTNNB1
selected(oncogenes[oncogenes$gene_name == "CTNNB1", "wmis_cv"],oncogenes[oncogenes$gene_name == "CTNNB1", "n_mis"])

##Para NFE2L2
selected(oncogenes[oncogenes$gene_name == "NFE2L2", "wmis_cv"],oncogenes[oncogenes$gene_name == "NFE2L2", "n_mis"])
  
## Para TP53
selected(data_genes[data_genes$gene_name == "TP53", "wmis_cv"],data_genes[data_genes$gene_name == "TP53", "n_mis"])
selected(data_genes[data_genes$gene_name == "TP53", "wnon_cv"],data_genes[data_genes$gene_name == "TP53", "n_non"])
selected(data_genes[data_genes$gene_name == "TP53", "wind_cv"],data_genes[data_genes$gene_name == "TP53", "n_ind"])

## Para AXIN1
selected(data_genes[data_genes$gene_name == "AXIN1", "wnon_cv"],data_genes[data_genes$gene_name == "AXIN1", "n_non"])
selected(data_genes[data_genes$gene_name == "AXIN1", "wind_cv"],data_genes[data_genes$gene_name == "AXIN1", "n_ind"])

## Para ARID1A
selected(data_genes[data_genes$gene_name == "ARID1A", "wnon_cv"],data_genes[data_genes$gene_name == "ARID1A", "n_non"])
selected(data_genes[data_genes$gene_name == "ARID1A", "wind_cv"],data_genes[data_genes$gene_name == "ARID1A", "n_ind"])

## Para KEAP1
selected(data_genes[data_genes$gene_name == "KEAP1", "wmis_cv"],data_genes[data_genes$gene_name == "KEAP1", "n_mis"])

## Para BAP1
selected(data_genes[data_genes$gene_name == "BAP1", "wnon_cv"],data_genes[data_genes$gene_name == "BAP1", "n_non"])

#Para ACVR2A
selected(data_genes[data_genes$gene_name == "ACVR2A", "wnon_cv"],data_genes[data_genes$gene_name == "ACVR2A", "n_non"])

# Para ARID2
selected(data_genes[data_genes$gene_name == "ARID2", "wnon_cv"],data_genes[data_genes$gene_name == "ARID2", "n_non"])

## Para TSC2
selected(data_genes[data_genes$gene_name == "TSC2", "wnon_cv"],data_genes[data_genes$gene_name == "TSC2", "n_non"])
selected(data_genes[data_genes$gene_name == "TSC2", "wind_cv"],data_genes[data_genes$gene_name == "TSC2", "n_ind"])

## Para RB1
selected(data_genes[data_genes$gene_name == "RB1", "wind_cv"],data_genes[data_genes$gene_name == "RB1", "n_ind"])

##Para NLGN1
selected(data_genes[data_genes$gene_name == "NLGN1", "wmis_cv"],data_genes[data_genes$gene_name == "NLGN1", "n_mis"])

#Para CRIBP3
selected(data_genes[data_genes$gene_name == "CRIP3", "wmis_cv"],data_genes[data_genes$gene_name == "CRIP3", "n_mis"])

#Para KRTAP3-3
selected(data_genes[data_genes$gene_name == "KRTAP3-3", "wnon_cv"],data_genes[data_genes$gene_name == "KRTAP3-3", "n_non"])

#Para ALB 
selected(data_genes[data_genes$gene_name == "ALB", "wind_cv"],data_genes[data_genes$gene_name == "ALB", "n_ind"])


###Total mutations
total_mut <- data_genes$n_syn + data_genes$n_mis + data_genes$n_non + data_genes$n_spl + data_genes$n_ind
newdata <- data.frame(data_genes$gene_name, total_mut)

## Hotspots

dout$annotmuts$gene_and_aachange = paste(dout$annotmuts$gene, dout$annotmuts$aachange, dout$annotmuts$ntchange, dout$annotmuts$pos, dout$annotmuts$impact, sep=":")
sort(table(dout$annotmuts$gene_and_aachange), decreasing = T)[1:15]



### Running sitednds
data("cancergenes_cgc81", package = "dndscv")
dout_cancergenes = dndscv(input, outmats = T, gene_list = known_cancergenes, max_muts_per_gene_per_sample = 3,max_coding_muts_per_sample = 500)
sd = sitednds(dout_cancergenes)
head(sd)
hotspots = sd$recursites[sd$recursites$qval < 0.05, ]
hotspots
hotspots[order(hotspots$qval), ]

### Ahora con codondnds
data("refcds_hg19", package = "dndscv")
RefCDS_codon = buildcodon(RefCDS)
codon_dnds = codondnds(dout_cancergenes, RefCDS_codon,
                       theta_option = "conservative",
                       min_recurr = 2)
codon_dnds$recurcodons[which(codon_dnds$recurcodons$qval < 0.1), ]

#############Oncoplots with maftools

library(maftools)
cancer = read.maf(maf = "~/LCG/2ndyear/3rd semester/Genomica funcional/data_mutations_maf.txt")
oncoplot(maf= cancer)
oncoplot(maf=cancer, top=30)
oncoplot(cancer, genesToIgnore = c("TTN"), top= 30)

genes_sig = c("TP53", "CTNNB1","ALB","AXIN1","ARID1A","KEAP1","TSC2","BAP1",
              "RB1","ACVR2A","NLGN1","CRIP3","KRTAP3-3","KRTAP22-1", "NFE2L2",
              "ARID2","IRX1")
oncoplot(maf = cancer,
         genes = genes_sig,
         sortByMutation = TRUE,
         removeNonMutated = TRUE,
         showTumorSampleBarcodes = FALSE)

oncoplot(cancer, genes = c("TP53","CTNNB1"))
oncoplot(cancer, genes = c("TP53","RB1"))
oncoplot(cancer, genes = c("KEAP1","NFE2L2"))
oncoplot(cancer, genes = c("CTNNB1","AXIN1"))


oncoplot(cancer, genes = genes_sig, draw_titv= TRUE)
oncoplot(maf=cancer, pathways= 'auto')
plotmafSummary(maf=cancer, addStat = 'median', dashboard = TRUE)
hepatocancer = titv(maf = cancer, plot=FALSE)
plotTiTv(res=hepatocancer)

### Lolipops
#TP53
lollipopPlot(maf = cancer,gene = "TP53",AACol = "HGVSp_Short",refSeqID = "NM_000546",showMutationRate = TRUE,labelPos = c(249, 193)) 
#CTNNB1
lollipopPlot(maf = cancer,gene = "CTNNB1",AACol = "HGVSp_Short",showMutationRate = TRUE,labelPos = c(45, 32, 36, 41, 33, 45))

