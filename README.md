# Probabilistic Mapping Extension

The **Probabilistic Mapping Extension** is a Python package designed to facilitate the creation of probabilistic maps for voxel-based datasets. It enables efficient processing of T-type and ME-type data by generating MET-type files, using probabilistic models to map relationships between input types.

---

## Features

- Probabilistic mapping of large datasets with ease.
- Handles `.nrrd` files for voxel density processing.
- Supports flexible configurations for mapping workflows.

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

## Usage Examples

Here are examples of how to use the **Probabilistic Mapping Extension**:

### Example 1: Basic Workflow

```python
from probabilistic_mapping import process_map

# Define input parameters
input_ttype = "path/to/ttype_file.nrrd"
input_pmap = "path/to/probability_map.csv"
output_met_type = "path/to/output_met_type.nrrd"

# Generate MET-type file
process_map(ttype=input_ttype, pmap=input_pmap, output=output_met_type)
print("Mapping completed!")
```

### Example 2: Advanced Configurations

```python
from probabilistic_mapping import process_map

# Custom configurations
process_map(
    ttype="path/to/ttype_file.nrrd",
    pmap="path/to/probability_map.csv",
    output="output_file.nrrd",
    threshold=0.8,  # Only process probabilities > 0.8
    verbose=True    # Enable detailed logs
)
```

---

## Contributing

Contributions are welcome! If you'd like to contribute:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
   ```bash
   git checkout -b my-feature-branch
   ```
3. Make your changes and commit them.
   ```bash
   git commit -m "Add my feature"
   ```
4. Push to your branch.
   ```bash
   git push origin my-feature-branch
   ```
5. Open a pull request on GitHub.

Please follow the guidelines in the [CONTRIBUTING.md](CONTRIBUTING.md) file.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

---

## Acknowledgments

Special thanks to all contributors and collaborators who helped bring this project to life!

# probabilistic_mapping_extention
