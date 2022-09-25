# set working directory to analysis data
setwd("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/analysis_data")

# download and load packages
install.packages("survival")
install.packages("survminer")
install.packages("ggplot2")
library(BiocManager)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(ggplot2)

# read in clinical dataset from clinical data tutorial
clinical <- read.csv("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/analysis_data/brca_clinical_data.csv")

# prepare drug and radiation data
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")

clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

# interesting variable from the clinical dataset: histological_type
sum(is.na(clinical$histological_type)) # 0 NAs

# interesting variable from the clinical.rad df: radiation_dosage
sum(is.na(clinical$lymph_node_examined_count)) # 139 NAs

# Coding Activity #1: basic boxplot comparing distribution of lymph nodes examined in different histological subtypes
# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/week5_hw/boxplot_histological_type_ln_examined_09262022.jpeg")
ggboxplot(clinical, x = "histological_type", y = "lymph_node_examined_count") +
  theme(axis.text.x = element_text(angle = 90))

# SURVIVAL ANALYSIS
# much of code for creating and fitting survival objects taken from Clinical_Data_Tutorial.Rmd

# defining variables
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)

clinical$death_event <- ifelse(clinical$vital_status == "Dead", TRUE, FALSE)

# survival analysis for the histological type variable
# Initialize a 'survival' object, which contains the data we need.
surv_object_histology <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
histology_fit <- surv_fit( surv_object_histology ~ clinical$histological_type,
                     data = clinical )

survplot_histology = ggsurvplot(histology_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/week5_hw/KM_histological_subtype_09262022.jpeg")
KM_plot_histology = survplot_histology$plot + 
  theme_bw() + # visual
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_histology

# survival analysis for the lymph nodes examined variable

# examine the rough distribution of the variable
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00    3.00    9.00   10.45   15.00   44.00     139
summary(clinical$lymph_node_examined_count)

# stratify patients by the # of lymph nodes they had examined using nested ifelse structure
# categories: 0-9, 10-19, 20-29, 30-39, 40+
clinical$ln_category <- ifelse(clinical$lymph_node_examined_count < 10, "0-9", ifelse(clinical$lymph_node_examined_count < 20, "10-19", ifelse(clinical$lymph_node_examined_count < 30, "20-29", ifelse(clinical$lymph_node_examined_count < 40, "30-39", "40+"))))

surv_object_ln <- Surv(time = clinical$survival_time,
                              event = clinical$death_event)

# Create a fit object
ln_fit <- surv_fit( surv_object_ln ~ clinical$ln_category, data = clinical )

survplot_ln = ggsurvplot(ln_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

# jpeg("/Users/brandonye/Desktop/fa22/QBIO 490/fa22_qbio_490_brandon/week5_hw/KM_ln_09262022.jpeg")
KM_plot_ln = survplot_ln$plot + 
  theme_bw() + # visuals
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_ln


