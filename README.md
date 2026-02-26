# Probabilistic Mapping Extension

The **Probabilistic Mapping Extension** is a Python package designed to facilitate the creation and analysis of probabilistic maps for voxel-based datasets. It provides tools to download experimental data, compute probabilistic maps, and perform optional density transcriptions and visualizations.

---

## Features

- Download experimental data for mapping workflows.
- Compute probabilistic maps in a structured sequence of steps.
- Optional transcription and visualization of voxel densities.
- Override specific regions for custom analyses.

---

## Installation

To install and use this package, follow these steps:

1. **Clone the Repository**  
   Clone the repository to your local machine:
   ```bash
   git clone https://github.com/YannRoussel/probabilistic_mapping_extention.git
   cd probabilistic_mapping_extention
   ```

2. **Install Dependencies**  
   Install the required dependencies with `pip`:
   ```bash
   pip install -r requirements.txt
   ```

3. **Set Up the Package**  
   Install the package locally using:
   ```bash
   pip install .
   ```

---

## Usage Guide

### Step 1: Download Experimental Data

Navigate to the `experimental_data` directory and run the scripts provided to download the necessary experimental data:
```bash
cd experimental_data
# Run the appropriate script(s)
```

### Step 2: Compute the Probabilistic Map

Sequentially execute the scripts in the `probability_map` directory:

1. Align T-types:
   ```bash
   python t_type_alignement.py
   ```

2. Compute the probabilistic ME-T map:
   ```bash
   python compute_p_me_t.py
   ```

3. Extend the ME-T map:
   ```bash
   python extend_p_me_t.py
   ```

### Step 3 (Optional): Transcribe and Visualize Densities

1. **Transcribe T-types to ME-types**  
   Use the `dens_app/transcribe_t_types_to_me_types.ipynb` Jupyter notebook to convert T-type densities (from another package) into ME-type densities.

2. **Visualize Results**  
   Visualize the results using the following notebooks:
   - `dens_app/Validation_compute_m_type_gene_expressions.ipynb`
   - `dens_app/visualisation_multiproc.ipynb`

### Step 4 (Optional): Override Specific Regions

Override certain ME-type densities using the `special_region_override/densities_override.py` script:
```bash
python special_region_override/densities_override.py
```

### Step 5 (Optional): Extension Folder

In this folder, you can find an extended analysis of the m/e/me cell types.

1. In the `data` folder, we generated cell type–by–brain area matrices for easier access to the densities.
   - All files are `.csv` files. We show how to access them in the notebook tutorials.
   - e-, m-, and me-types are kept separate, as suggested by the file names.
   - Brain region naming follows the Allen Institute nomenclature. You can include unassigned regions by selecting the file that contains this in its file name.
   - `all_celltypes` means everything is combined into a single large file.

2. In the `etype` folder, you can find generated figure insets for all (ephys) “electric types,” shown as cross-sections of generated axial or transversal example brain slices. These are located in:
   - `etype/fig_me_ax` (axial)
   - `etype/fig_me_tr` (sagittal)

   The notebooks demonstrate how to generate and concatenate these figures.

3. In the `figures` folder, you can find all m/e/me-type insets concatenated into composite images. The notebooks show how to generate these, including the brain slices as 2D arrays in `.npy` format (not included here).

4. In the `laminar_validation` folder, you can find code to generate cell type compositions plotted across different brain regions, along with the resulting images (in `laminar_validation/results`). As an example, we plot excitatory and inhibitory m/e/me-types in the AUD, VISp, and SSs layers of the CTX.

5. In the `me-types` folder, you can find generated figure insets for all “morpho-electric types,” shown as cross-sections of generated axial or transversal example brain slices. These are located in:
   - `me-types/fig_me_ax` (axial)
   - `me-types/fig_me_tr` (sagittal)

   The notebooks demonstrate how to generate and concatenate these figures.

6. In the `mtypes` folder, you can find generated figure insets for all “morphological types,” shown as cross-sections of generated axial or transversal example brain slices. These are located in:
   - `m-types/fig_me_ax` (axial)
   - `m-types/fig_me_tr` (sagittal)

   The notebooks demonstrate how to generate and concatenate these figures.

7. The `notebooks` folder contains Jupyter notebook tutorials:
   - `compare_your_data.ipynb`: Visualize any cell type in any brain slice together with your own data. This is the only notebook that uses `ipywidgets`.
   - `create_all_celltype_df.ipynb`: Intermediate step showing how to concatenate cell types into `all_celltypes_composition_with_unassigned_regions.csv` in the `data` folder.
   - `e_type_slices_for_figure.ipynb`: Generate `.npy` files (2D NumPy arrays) of brain slices with e-type densities.
   - `laminarity_structure_metypes.ipynb`: Test code demonstrating laminar organization of me-types for validation. This notebook helps clarify the corresponding Python implementation.
   - `laminarity_structure_mtypes.ipynb`: Test code demonstrating laminar organization of m-types for validation.
   - `m_type_slices_for_figure.ipynb`: Generate `.npy` files (2D NumPy arrays) of brain slices with m-type densities.
   - `me_type_slices_for_figure.ipynb`: Generate `.npy` files (2D NumPy arrays) of brain slices with me-type densities.
   - `transcribing_t_types.ipynb`: Simplified version of `densities_app/transcribe_t_types_to_me_types.ipynb`, focusing only on the necessary steps and including instructions on saving `.csv` files in the `data` folder.
   - `visualisation_multiproc.ipynb`: Simplified version of `densities_app/visualisation_multiproc.ipynb`, focusing only on the necessary steps.
   - `visualize_celltypes_supp_figures.ipynb`: Notebook showing how to generate selected supplementary figures for the white paper. The resulting composite images are available in the `figures` folder.

8. `annotation_25_2022_CCFv3a.nrrd`: All notebooks use this Common Coordinate Framework (from the Blue Brain Cell Atlas), which is a modified version of the Allen Institute’s CCFv3.


---

## Contributing

Contributions are welcome! If you'd like to contribute:

1. Fork the repository.
2. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b my-feature-branch
   ```
3. Make your changes and commit them:
   ```bash
   git commit -m "Add my feature"
   ```
4. Push to your branch:
   ```bash
   git push origin my-feature-branch
   ```
5. Open a pull request on GitHub.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

## Acknowledgments

Special thanks to all contributors and collaborators who helped bring this project to life!

