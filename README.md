# PIEZO1 puncta localization analysis pipeline

## Git Clone repository to your host system

```
git clone --recurse-submodules https://github.com/abcucberkeley/piezo1_analysis_pipeline.git

# To later update the repository:
git pull --recurse-submodules
```

## Usage

* The codes have been tested with MATLAB R2022b-R2023a for Linux (Ubuntu).

* The `main_pipeline/` directory contains three scripts:
    *  `puncta_localization_analysis_2d.m`: This script generates svg files of density scatter and frequency distribution plots for the distances of puncta localizations to Lumen Edge and Outer Edge Masks, along with csv files containing cdf values of the distances.
    *  `puncta_tracking_analysis_2d.m`: This script generates svg files of density scatter and frequency distribution plots for the distances of the mean position of puncta tracks to Lumen Edge and Outer Edge Masks, along with csv files containing cdf values of the distances.
    *  `filtered_synthetic_3d_volume_generation.m`: This script generates 3D synthetic volumes for PIEZO1 puncta detections and autofluorescence blobs obtained using `LLSM5DTools/`, and the local density map of the PIEZO1 detections. The 3D volumetric data is then used for visualization using **Imaris x64: 10.0.0** and **Amira 3D 2021.1**.

* For each of the scripts, the corresponding results are stored inside the root directory.
    
* The scripts for the 2D analysis require the root directory to contain separate subdirectories for each individual video. Each subdirectory should contain a "results" folder, which stores the csv files for puncta localization and tracking obtained from **ThunderSTORM**, and an "actin_mask" folder, which stores the tif files for both the lumen and outer edge masks.

* The density scatter and frequency distribution plots shown in the paper [PIEZO1-HaloTag hiPSCs: Bridging Molecular, Cellular and Tissue Imaging](https://doi.org/10.1101/2023.12.22.573117) can be replicated by downloading the dataset from https://doi.org/10.5061/dryad.w6m905qwm, and running the scripts in the `main_pipeline/` directory, after setting the correct path for rootdir.

* All functional dependencies are located in the `util_functions/` directory and `LLSM5DTools/` submodule.

* The `nuclei_denoising_model/` directory contains a deep learning model trained using **CARE** for denoising the nuclei channel in the 3D volumetric data.
