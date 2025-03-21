# Whole-brain 3D maps of cellular abundance for six canonical cell-types in humans

<p align="center">
<img src="img.png" width="600" height="300">
</p>

High-resolution imputed maps of spatial cell-type distributions can be used to study variations in cell-type density in relation to macroscale phenotypes. These maps can be utilized in studies focused on both the structure and function of the brain in health and disease, specifically within the context of neurological conditions including neurodegenerative, neurodevelopmental, and psychiatric disorders.


(c) Neuroinformatics for Personalized Medicine Lab (Neuro-PM): https://www.neuropm-lab.com/

Authors: Veronika Pak, Quadri Adewale, Professor Yasser Iturria Medina

DOI: 10.1101/2023.06.08.544227

Contact information: veronika.pak@mail.mcgill.ca

## Maps Description
- six brain cell-types: ast - astrocytes, end - endothelial cells, mic - microglia, neu - neurons, oli - oligodendrocytes, opc - oligodendrocyte precursor cells
- NIfTI format
- 1.5 mm resolution
- ~460,600 gray matter voxels across the whole brain
- registered to the Montreal Neurological Institute (MNI) brain space
- normalized and scaled relative to gray matter volume distribution
- reconstructed from post-mortem bulk gene expression data from the Allen Human Brain Atlas
- for methods, see [publication](https://doi.org/10.1101/2023.06.08.544227)

## regionalCorr.m
The MATLAB code allows you to extract regional values, perform statistical analysis, and compute correlations between the regional cell type proportions and the neuroimaging phenotype across the atlas regions.

### Features

- **Flexible Correlation Methods:**  
  Choose `'Spearman'` (default) or `'Pearson'` as the correlation method.

- **Automatic Data Cleaning:**  
  Automatically removes regions with missing (0 or NaN) data and reports which rows were removed.

### Requirements

- **MATLAB:**  
  Implemented in MATLAB.

- **SPM:**  
  Uses `spm_vol` for converting voxel coordinates to MNI space. Ensure that SPM is installed and added to your MATLAB path.

- **NIfTI Support:**  
  Relies on MATLAB's `niftiread` function for reading NIfTI files.

### Usage Example

```matlab
% Define the file names
atlasFile = 'DTKatlas.nii'; % Change to your parcellation of interest
xFiles = {'ast.nii', 'end.nii', 'mic.nii', 'neu.nii', 'oli.nii', 'opc.nii'};
yFile = 'LOAD_atrophy.nii'; % Change to your phenotype of interest

% Compute correlations using Spearman correlation (change to 'Pearson if needed)
[rho, pval, regionList, xRegional, yRegional] = regionalCorr(atlasFile, xFiles, yFile, 'Spearman');

% Display results
disp('Correlation coefficients for each cell type:');
disp(rho);
disp('P-values:');
disp(pval);
```

## Citation
``` Veronika Pak, Quadri Adewale, Danilo Bzdok, Mahsa Dadar, Yashar Zeighami, Yasser Iturria-Medina. Distinctive Whole-brain Cell-Types Patterns Strongly Predict Tissue Damage in Eleven Neurodegenerative Disorders. eLife. https://elifesciences.org/articles/89368 ```
