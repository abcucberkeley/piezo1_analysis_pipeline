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
    *  `puncta_localization_analysis_2d.m`: This script generates plots for the distances of puncta localizations to Lumen Edge and Outer Edge Masks.
    *  `puncta_tracking_analysis_2d.m`: This script generates plots for the distances of the mean position of puncta tracks to Lumen Edge and Outer Edge Masks.
    *  `filtered_synthetic_3d_volume_generation.m`: This script generates 3D synthetic volumes for PIEZO1 puncta detections and autofluorescence blobs obtained using `LLSM5DTools/`, and the local density map of the PIEZO1 detections.
    
* The scripts for the 2D analysis require the root directory to contain separate subdirectories for each individual video, which in turn contains the csv files for puncta localization and tracking obtained from **ThunderSTORM**, and the tif files for both the masks. 

* All functional dependencies are located in the `util_functions/` directory and `LLSM5DTools/` submodule.

* The `nuclei_denoising_model/` directory contains a deep learning model trained using **CARE** for denoising the nuclei channel.
