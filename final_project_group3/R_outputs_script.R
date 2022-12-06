---
title: "R Notebook"
output: html_notebook
---

if (!require(TAF)) {install.packages("TAF")}
library(TAF)
setwd("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon")
mkdir("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/final_project_group3")
mkdir("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/final_project_group3/outputs")
setwd("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/final_project_group3/outputs")

install.packages("TCGAbiolinks")
install.packages("maftools")
install.packages("BiocManager")
install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")
#install.packages("vioplot")
#install.packages("aplpack")

library("TCGAbiolinks")
library("maftools")
library("BiocManager")
library("survival")
library("survminer")
library("ggplot2")
#library("vioplot")
#library("aplpack")

clinical_query <- GDCquery(project = "TCGA-LUAD", data.category = "Clinical", file.type = "xml")
GDCdownload(clinical_query)

clinical <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient")

clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#filters N/A values from dataframe
clinical_filtered = clinical[!is.na(clinical$year_of_initial_pathologic_diagnosis), ]

#creates smoking_status column in dataframe based on whether there is data on year of first smoke
clinical_filtered$smoking_status =
  ifelse(is.na(clinical_filtered$year_of_tobacco_smoking_onset), FALSE, TRUE)

#graphs boxplot of age at LA diagnosis based on smoking status
jpeg("clinical_boxplot.jpg")
par(mar=c(4,4,4,4))
boxplot(clinical_filtered$age_at_initial_pathologic_diagnosis ~ clinical_filtered$smoking_status, xlab = "Smoking Status", ylab = "Age of Diagnosis", main = "Age of Lung Adenocarcinoma Diagnosis\nBased on Smoking Status")
dev.off()

#graphs histogram of difference in years from first smoke to LA diagnosis
jpeg("clinical_histogram.jpg")
par(mar=c(4,4,4,4))
hist(clinical_filtered$year_of_initial_pathologic_diagnosis -
       clinical_filtered$year_of_tobacco_smoking_onset, xlab = "Years", ylab = "Frequency", main = "Years from Start of Smoking to\nInitial Lung Adenocarcinoma Diagnosis")
dev.off()

#download all maf files and queried in
maf_query <- GDCquery(
  project = "TCGA-LUAD", data.category = "Simple Nucleotide Variation", access = "open",
  data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf <- GDCprepare(maf_query) # as long as it runs, ignore any errors

colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_object <- read.maf(maf = maf, clinicalData = clinical, isTCGA = TRUE)

#creating stratification to determine smoking status
clinical$smoking_status = ifelse(is.na(clinical$year_of_tobacco_smoking_onset), FALSE, TRUE)

#dataframe just with smoking patients
clinical_smoking <- clinical[clinical$smoking_status == TRUE,]

#stratfication to determine survival time and death status, copied onto smoking population
clinical_smoking$survivaltime = ifelse(is.na(clinical_smoking$days_to_death), clinical_smoking$days_to_last_followup, clinical_smoking$days_to_death)
clinical_smoking$deathstatus = ifelse(clinical_smoking$vital_status == "Dead", TRUE, FALSE)

#Kaplan Meier Survival Analysis of Gender within Smoking Patients
survival_object <- Surv(time = clinical_smoking$survivaltime, event = clinical_smoking$deathstatus)

gender_sort <- surv_fit(survival_object ~ clinical_smoking$gender, data = clinical_smoking)
jpeg("survival_sex.jpg")
ggsurvplot(gender_sort, data = clinical_smoking, surv.median.line = "hv", legend.title = "Gender", legend.labs = c("Men", "Women"), pval = TRUE, conf.int = TRUE, risk.table = TRUE, tables.height = 0.2, tables.theme = clean_theme(), ggtheme = theme_gray())
dev.off()

#stratfication to determine survival time and death status, copied onto smoking population
clinical$survivaltime = ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)
clinical$deathstatus = ifelse(clinical$vital_status == "Dead", TRUE, FALSE)

#Kaplan Meier Survival Analysis of Smoking Patients
survival_object <- Surv(time = clinical$survivaltime, event = clinical$deathstatus)
smoking_sort <- surv_fit(survival_object ~ clinical$smoking_status, data = clinical)

jpeg("survival_smoking.jpg")
ggsurvplot(smoking_sort, data = clinical, surv.median.line = "hv", legend.title = "Smoking Status", legend.labs = c("TRUE", "FALSE"), pval = TRUE, conf.int = TRUE, risk.table = TRUE, tables.height = 0.2, tables.theme = clean_theme(), ggtheme = theme_gray())
dev.off()

#read in updated maf object
maf_object <- read.maf(maf = maf, clinicalData = clinical, isTCGA = TRUE)

#draw oncoplot with smoking patients
jpeg("top_10_oncoplot.jpg")
oncoplot(maf = maf_object, top = 10, clinicalFeatures = "smoking_status")
dev.off()

#cleaning maf data to differentiate between smoking positive and smoking negative patient barcodes
smoking_positive_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$smoking_status == TRUE]

smoking_positive_maf <- subsetMaf(maf= maf_object, tsb = smoking_positive_barcodes)

smoking_negative_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$smoking_status == FALSE]

smoking_negative_maf <- subsetMaf(maf = maf_object, tsb = smoking_negative_barcodes)

#draw cooncoplot comparing smokers and nonsmoker gene mutation frequencies
jpeg("cooncoplot_smoking.jpg")
coOncoplot(m1 = smoking_positive_maf, m2 = smoking_negative_maf, m1Name = "Smoker", m2Name = "Non-Smoker")
dev.off()

#map gene locus of the gene with the greatest differential mutation rate within the top 5 most mutated genes, which in this case is TTN
jpeg("lollipop_smoking_TTN.jpg")
lollipopPlot2(m1 = smoking_positive_maf, m2 = smoking_negative_maf, m1_name = "Smoker", m2_name = "Non-Smoker", gene = "TTN")
dev.off()

#secondary analysis of sex within smoking positive: sorting patient sex by barcode
male_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$gender == "MALE"]

male_maf <- subsetMaf(maf= smoking_positive_maf, tsb = male_barcodes)

female_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$gender == "FEMALE"]

female_maf <- subsetMaf(maf = smoking_positive_maf, tsb = female_barcodes)

#draw cooncoplot to compare rates of mutation for top 5 genes between sex stratification of smoking cohort
jpeg("cooncoplot_sex.jpg")
coOncoplot(m1 = male_maf, m2 = female_maf, m1Name = "Male", m2Name = "Female")
dev.off()

#lollipop plot to map gene locus of CSMD3, which had the greatest differential mutation rate within the sex stratification of the smoking cohort
jpeg("lollipop_sex_CSMD3.jpg")
lollipopPlot2(m1 = male_maf, m2 = female_maf, m1_name = "Male", m2_name = "Female", gene = "CSMD3")
dev.off()

#download and query of RNA data
rna_query <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)

#masks on RNA data to create rna counts, rna genes, and rna clinical data frames for analysis
age_mask <-  is.na(rna_se@colData$age_at_index)
rna_counts <- rna_se@assays@data$unstranded[, !age_mask]

rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)

rna_clinical <-  rna_se@colData[!age_mask, ]

#masks for each of the highest frquency gene mutations within the TCGA-LUAD database
TP53_mask <- ifelse(rna_genes$gene_name == "TP53", TRUE, FALSE)
TP53_counts <- as.numeric(rna_counts[TP53_mask, ])

KRAS_mask <- ifelse(rna_genes$gene_name == "KRAS", TRUE, FALSE)
KRAS_counts <- as.numeric(rna_counts[KRAS_mask, ])

STK11_mask <- ifelse(rna_genes$gene_name == "STK11", TRUE, FALSE)
STK11_counts <- as.numeric(rna_counts[STK11_mask, ])

EGFR_mask <- ifelse(rna_genes$gene_name == "EGFR", TRUE, FALSE)
EGFR_counts <- as.numeric(rna_counts[EGFR_mask, ])

LRP1B_mask <- ifelse(rna_genes$gene_name == "LRP1B", TRUE, FALSE)
LRP1B_counts <- as.numeric(rna_counts[LRP1B_mask, ])

CSMD3_mask <- ifelse(rna_genes$gene_name == "CSMD3", TRUE, FALSE)
CSMD3_counts <- as.numeric(rna_counts[CSMD3_mask, ])

TTN_mask <- ifelse(rna_genes$gene_name == "TTN", TRUE, FALSE)
TTN_counts <- as.numeric(rna_counts[TTN_mask, ])

MUC16_mask <- ifelse(rna_genes$gene_name == "MUC16", TRUE, FALSE)
MUC16_counts <- as.numeric(rna_counts[MUC16_mask, ])

RYR2_mask <- ifelse(rna_genes$gene_name == "RYR2", TRUE, FALSE)
RYR2_counts <- as.numeric(rna_counts[RYR2_mask, ])

#sorting data and survival data to determine survival and smoking status within RNA data
geneABCD_counts <- data.frame(TP53_counts, KRAS_counts, STK11_counts, EGFR_counts, LRP1B_counts, TTN_counts, MUC16_counts, CSMD3_counts, RYR2_counts)
colnames(geneABCD_counts) <- c("TP53", "KRAS", "STK11", "EGFR", "LRP1B", "TTN", "MUC16", "CSMD3", "RYR2")

five_yr_death <- ifelse(rna_clinical$days_to_death == "NA", NA, ifelse(rna_clinical$days_to_death > 365.25*5, TRUE, FALSE))
five_yr_death_and_followup <- ifelse(is.na(five_yr_death), ifelse(rna_clinical$days_to_last_follow_up > 365.25*5, TRUE, FALSE), five_yr_death)
rna_clinical$five_year_surv <- five_yr_death_and_followup

smoking_status <- ifelse(rna_clinical$years_smoked == "NA", FALSE, TRUE)

rna_clinical$smoking_status <- ifelse(is.na(rna_clinical$years_smoked), FALSE, TRUE)

# these lines set up a clustering color scheme for our plot, chose smoking status from rna_clinical to cluster based on 
cols <- character(nrow(rna_clinical)) 
cols[rna_clinical$smoking_status == TRUE] <- "blue" 
cols[rna_clinical$smoking_status == FALSE] <- "red"

# draftsman plot code
jpeg("draftsman.jpg")
pairs(geneABCD_counts, col = cols, lower.panel=NULL)
dev.off()
