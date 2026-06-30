# Tensor-Based HO-MSSA for Biomedical Signal Analysis

This repository contains the MATLAB implementation, helper functions, and datasets for the paper: **"Tensor-based higher-order multivariate singular spectrum analysis and applications to multichannel biomedical signal analysis."** The code implements Higher-Order Multivariate Singular Spectrum Analysis (HO-MSSA) via tensor decomposition to process and reconstruct multichannel time series signals, specifically for biomedical applications like ECG and EEG.

## 📄 Reference Paper
If you use this code or data in your research, please cite the following paper:
> Thanh Trung Le, Karim Abed-Meraim, Nguyen Linh Trung, Philippe Ravier, Olivier Buttelli, Ales Holobar. *"Tensor-based higher-order multivariate singular spectrum analysis and applications to multichannel biomedical signal analysis."* Signal Processing, Volume 238 (2026), 110113.

---

## 🛠️ Installation and Setup

To run this code, simply download all the files in this repository and place them together in a single folder on your machine or MATLAB Online environment. Ensure that this folder is set as your Current Folder in MATLAB.

---

## 📂 Repository Structure

The files in this repository are divided into three main categories: Data, Main Scripts, and Helper Functions.

### 1. Data Files
* `FOETAL_ECG.mat`: Electrocardiogram (ECG) dataset.
* `Contaminated_Data.mat`: Electroencephalogram (EEG) dataset.

### 2. Main MATLAB Scripts
These are the primary executable scripts for running the analyses:

* **`HOMSSA.m`**
    * *Description:* [Add your details here about what this script does, inputs/outputs, etc.]
* **`batch_evaluate_table2.m`**
    * *Description:* [Add your details here about what this script does, inputs/outputs, etc.]
* **`calculateK.m`**
    * *Description:* [Add your details here about what this script does, inputs/outputs, etc.]

### 3. Helper Functions
These are underlying functions called by the main scripts to perform tensor operations, clustering, and signal reconstruction:
* `tblock_diagonal_averaging.m`
* `tbuild_trajectory_tensor.m`
* `tcluster_components.m`
* `tcompute_dynamic_tubal_svd.m`
* `tcompute_tubal_metrics.m`
* `tevaluate_mae_frequency.m`
* `tplot_figure_13_single.m`
* `treconstruct_signals.m`
* `tsvd.m`
* `tubalrank.m`

---

## 🚀 Usage

[Add any general instructions here on how a user should begin. For example: "To reproduce the results from Table 2 in the paper, run `batch_evaluate_table2.m` directly in the command window."]