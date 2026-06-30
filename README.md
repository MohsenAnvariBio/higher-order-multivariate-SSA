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
    * *Description:* This is the main executable script that implements the 5-step Higher-Order Multivariate Singular Spectrum Analysis (HO-MSSA) algorithm. It processes the input data through Time Delay Embedding (TDE), Tensor SVD, Spectral Clustering, Signal Reconstruction, and finally generates waveform visualizations.
    * *Data Selection (EEG vs. ECG):* You can easily switch between processing the EEG or ECG datasets by commenting or uncommenting the respective `%% EEG DATA` or `%% ECG DATA` blocks at the very beginning of the script (Step 0). 
        * **Note for EEG:** When running the EEG data block, the script will prompt you in the Command Window to enter a patient ID (1 to 27) to dynamically load and process that specific participant's dataset.
    * *Algorithm Flow:*
        1. **Construct Trajectory Tensor:** Embeds the multichannel time series into a tensor.
        2. **Tubal SVD:** Decomposes the tensor dynamically using a 99.9% variance threshold to find the optimal tubal rank (reproducing Figure 13 from the paper).
        3. **Grouping:** Clusters the components into distinct groups (e.g., source signals vs. artifacts) using spectral clustering.
        4. **Reconstruction:** Uses block diagonal averaging (De-TDE) to recover the time-domain signals.
        5. **Visualization:** Automatically plots the original mixture alongside the separated clusters for each channel.

* **`calculate_mae.m`**
    * *Description:* [Add your details here about what this script does, inputs/outputs, etc.]
* **`calculate_K.m`**
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