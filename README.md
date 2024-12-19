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

