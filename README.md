# PIEZO1 puncta localization analysis pipeline

## Git Clone repository to your host system

```
git clone --recurse-submodules https://github.com/abcucberkeley/piezo1_analysis_pipeline.git

# To later update the repository:
git pull --recurse-submodules
```

## Hardware

The [demos](https://github.com/abcucberkeley/piezo1_analysis_pipeline/tree/main/main_pipeline) have been tested on systems with Intel Core i7 processor and 16 GB of RAM, as well as Intel Core i9 processor and 32 GB of RAM. We recommend at least 16 GB of RAM, and GPU is not required.

## Usage

The codes have been tested with **MATLAB R2022b-R2023b** for Linux (Ubuntu), Windows 10 and MacOS Ventura. The **Image Processing Toolbox** and **Statistics and Machine Learning Toolbox** for MATLAB are required for running the codes.

The scripts for the 2D analysis require the root directory to contain separate subdirectories for each individual video. Each subdirectory should contain a `results/` folder, which stores the csv files for puncta localization and tracking obtained from **ThunderSTORM**, and an `actin_mask/` folder, which stores the tif files for both the lumen and outer edge masks.

The density scatter and frequency distribution plots shown in the paper [PIEZO1-HaloTag hiPSCs: Bridging Molecular, Cellular and Tissue Imaging](https://doi.org/10.1101/2023.12.22.573117) can be replicated by following the steps: 
* Download the dataset from the [Dryad repository](https://doi.org/10.5061/dryad.w6m905qwm), and unzip it to a directory.
* Get the source code by either cloning the [GitHub repository](https://github.com/abcucberkeley/piezo1_analysis_pipeline.git), or downloading the ZIP file and unzipping the file to a directory.
* Launch MATLAB, navigate to the `main_pipeline/` directory, and run the codes after setting the correct path for `rootdir`. 
* The plot parameters such as **Position** might need to be adjusted depending on the display resolution.
* For each of the scripts, the corresponding results are stored inside the `results/` subdirectory within `rootdir`.

## Directory Structure

* The `main_pipeline/` directory contains three demo scripts:
    *  `puncta_localization_analysis_2d.m`: This script generates svg files of density scatter and frequency distribution plots for the distances of puncta localizations to Lumen Edge and Outer Edge Masks, along with csv files containing cdf values of the distances.
    *  `puncta_tracking_analysis_2d.m`: This script generates svg files of density scatter and frequency distribution plots for the distances of the mean position of puncta tracks to Lumen Edge and Outer Edge Masks, along with csv files containing cdf values of the distances.
    *  `filtered_synthetic_3d_volume_generation.m`: This script generates 3D synthetic volumes for PIEZO1 puncta detections and autofluorescence blobs obtained using `LLSM5DTools/`, and the local density map of the PIEZO1 detections. The 3D volumetric data is then used for visualization using **Imaris x64: 10.0.0** and **Amira 3D 2021.1**.

* All functional dependencies are located in the `util_functions/` directory and `LLSM5DTools/` submodule.

* The `nuclei_denoising_model/` directory contains a deep learning model trained using **CARE** for denoising the nuclei channel in the 3D volumetric data.
