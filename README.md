# Xenium Anndata Processing and Analysis Pipeline

This repository provides a step-by-step pipeline for processing Xenium data, merging datasets, and performing clustering using Banksy and Harmony.

---

## **Pipeline Overview**

### **Step 1: Generate Anndata for Each Xenium Sample**
Use the provided script to generate anndata objects from individual Xenium samples.

Run the following command:
```bash
bash run_making_anndata.sh
```
### **Step 2: Merge anndata and Stagger the Coordinates**
Merge the generated anndata objects and stagger their spatial coordinates. Optionally, remove unneeded genes by specifying a gene list.

Run the following command:
```bash
bash run_merge_anndata.sh
```
### **Step 3: Run Banksy and Harmony, Then Perform Clustering**
Run Banksy and Harmony for data integration, followed by clustering.

Run the following command:
```bash
bash run_Banksy.sh
```
