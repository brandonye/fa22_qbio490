# Mid-Semester Project
# QBIO 490 Fall 2022
# Brandon Ye

install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")
install.packages("maftools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.15")
library(BiocManager)
  
if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

library(maftools)
library(survival)
library(survminer)
library(ggplot2)

knitr::opts_knit$set(root.dir = normalizePath("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs")) 

# query for clinical data
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient")

# prepare (rename) column names of clinical for our maf query
colnames(clinical)[colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

# prepare clinical.drug and clinical.rad data
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

# query for MAF data
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(maf_query)

maf <- GDCprepare(maf_query) 

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

# the aim of this project is to determine how race affects survival in breast cancer patients, with a focus on differential mutation rates in TP53 and PIK3CA by race, and how these rate ultimately affect survival outcomes

# create a clinical survival time column for KM analysis
# if days_to_death is not NA, that is survival time
# if days_to_death is NA, days_to_last_followup becomes the best estimate
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)

# create a death event column; T if death recorded, F if otherwise
clinical$death_event <- ifelse(clinical$vital_status == "Dead", T, F)

# data cleaning for KM visualization
unique(clinical$race_list) # the categories that do exist look to be clean, however, there are many patients with no race information
clinical$race_list <- sub("^$", "OTHER", clinical$race_list) # replace patients who have no race information with race "OTHER"

# start with a simple countplot to visualize the distribution of the clinical dataset by race
# AMERICAN INDIAN OR ALASKA NATIVE                            ASIAN        BLACK OR AFRICAN AMERICAN 
#                                1                               62                              201 
#                            OTHER                            WHITE 
#                              109                              801 

# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/countplot_race_clinical.jpeg")
countplot_race_clinical <- ggplot(clinical, aes(x=race_list))
countplot_race_clinical + geom_bar()
# dev.off()

# plot a KM to depict differential survival by race
# create a survival object for race
surv_object_race <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
# Plot the surv_object against the race
race_fit <- surv_fit(surv_object_race ~ clinical$race_list,
                     data = clinical)

# plot formatting
survplot_race = ggsurvplot(race_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# finalize KM plot
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/KM_plot_race.jpeg")
KM_plot_race = survplot_race$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_race
# dev.off()

# merge clinical and clinical.drug
colnames(clinical)[colnames(clinical) == "Tumor_Sample_Barcode" ] <- "bcr_patient_barcode" # temp rever tumor sample barcode column name to merge
clinical_and_drug <- merge(clinical, clinical.drug, by = "bcr_patient_barcode", all.x=TRUE, all.y=TRUE)
colnames(clinical)[colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

unique(clinical_and_drug$therapy_types) # lots of "other" types which we should combine into one

# combine all the "Other" types... I do not know how to use regex in R
clinical_and_drug$therapy_types <- ifelse(clinical_and_drug$therapy_types == "Chemotherapy", "Chemotherapy", 
                                          ifelse(clinical_and_drug$therapy_types == "Hormone Therapy", "Hormone Therapy", 
                                                 ifelse(clinical_and_drug$therapy_types == "Vaccine", "Vaccine", 
                                                        ifelse(clinical_and_drug$therapy_types == "Immunotherapy", "Immunotherapy",
                                                               ifelse(clinical_and_drug$therapy_types == "Ancillary", "Ancillary", "Other")))))

# make a grouped countplot that depicts the counts for each therapy type by race
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/countplot_race_therapy.jpeg")
countplot_race_therapy <- ggplot(clinical_and_drug, aes(x=race_list, fill=therapy_types)) + geom_bar(position="dodge") # construct countplot
countplot_race_therapy + ylab("Count") + scale_fill_discrete(name="therapy types", labels = c("Chemotherapy", "Hormone Therapy", "Immunotherapy", "Vaccine", "Ancillary", "Other", "NA")) # legend and y-axis
# dev.off()

# TP53 and PIK3CA are the commonly mutated gene in BRCA patients, and are understood to possibly have antagonistic effects on survival
# Here, we move from clinical to MAF data and investigate how, if at all, race is associated with mutation rates in TP53 and PIK3CA

# a basic oncoplot for TP53 and PIK3CA
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/TP53_PIK3CA_oncoplot.jpeg")
genes <- c("TP53", "PIK3CA")
oncoplot(maf = maf_object,
         genes = genes)
# dev.off()

# we can use a CoOncoplot to visualize the differential distribution in these two genes by race
# to do so, we need to subset our maf into separate sections by race
# because CoOncoplots only allow two plots to displayed side by side, we will distinguish race by white (majority class) versus non-white (minority class)
# unique(maf_object@clinical.data$race_list) # WHITE, ASIAN, AMERICAN INDIAN OR ALASKA NATIVE, BLACK OR AFRICAN AMERICAN
white_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$race_list == "WHITE"]
other_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[maf_object@clinical.data$race_list != "WHITE"]

# subset the overall MAF with the divisions we just created
white_maf <- subsetMaf(maf = maf_object,
                       tsb = white_patient_barcodes)

other_maf <- subsetMaf(maf = maf_object,
                       tsb = other_patient_barcodes)

# plot using coOncoplot()
# whites have a lower mutation rate in TP53, and a higher mutation rate in PIK3CA
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/race_coOncoplot.jpeg")

coOncoplot(m1 = white_maf, m2 = other_maf, m1Name = "White", m2Name = "Other", genes = genes)

# dev.off()

# Co-lollipop plots are useful because they can depict mutations at the genome level
# is the differential mutation rate in TP53 and PIK3CA noticeable on the genome level?

# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/TP53_race_coLollipop.jpeg")
lollipopPlot2(m1 = white_maf, 
              m2 = other_maf, 
              m1_name = "White",
              m2_name = "Other",
              gene = "TP53") 
# dev.off()

# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/PIK3CA_race_coLollipop.jpeg")
lollipopPlot2(m1 = white_maf, 
              m2 = other_maf, 
              m1_name = "White",
              m2_name = "Other",
              gene = "PIK3CA") 
# dev.off()

# the differential mutation rate in TP53 and PIK3CA seem to be unnoticeable on the genome level
# the final analysis will involve determining whether or not mutations in each gene affects survival to a significant extent

# recreate the overall survival status and survival time columns for the mafSurvival plots
maf_object@clinical.data$overall_survival_status <- ifelse(maf_object@clinical.data$vital_status == "Dead", TRUE, FALSE)
maf_object@clinical.data$survival_time <- ifelse(is.na(maf_object@clinical.data$days_to_death), maf_object@clinical.data$days_to_last_followup, maf_object@clinical.data$days_to_death)

# Kaplan Meier mutated/non-mutated for TP53
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/TP53_mafSurvival.jpeg")
mafSurvival(maf = maf_object,
            genes = "TP53", 
            time = "survival_time", 
            Status = "overall_survival_status", 
            isTCGA = TRUE)
# dev.off()

# Kaplan Meier mutated/non-mutated for PIK3CA
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/midsemester_project_ye/outputs/PIK3CA_mafSurvival.jpeg")
mafSurvival(maf = maf_object,
            genes = "PIK3CA", 
            time = "survival_time", 
            Status = "overall_survival_status", 
            isTCGA = TRUE)
# dev.off()


