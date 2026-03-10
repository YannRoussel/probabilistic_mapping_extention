# SCALA:
### in 'probabilistic_mapping_extension/experimental_data/scala/' execute

## metadata
https://download.brainimagelibrary.org/3a/88/3a88a7687ab66069/m1_patchseq_meta_data.csv  

## Electrophysiology:
### main dataset:
dandi download DANDI:000008/0.211014.0809
### physiological dataset:
dandi download DANDI:000035/0.211014.0808

## Morphologies:
python morphologies/dl_scala_morphologies.py

## scRNAseq:
python rna/dl_rna_data_scala.py

# GOUWENS:
### in 'probabilistic_mapping_extension/experimental_data/gouwens/' execute

## metadata
./20200625_patchseq_metadata_mouse.csv
./20200711_patchseq_metadata_mouse(1).csv

## Electrophysiology
dandi download DANDI:000020/0.210913.1639

## Morphologies
python morphologies/dl_gouwens_morphologies.py

## scRNAseq
python rna/dl_rna_data_gouwens.py

# Yao et al. realigned
### obtain the t-types cluster 'median rna expression' and 'trimmed mean rna expression' by running the pipeline at https://github.com/BlueBrain/Molsys-transcriptomic-atlas

# canonical m_types (Kanari et al.)

./canonical m-types/test_cluster.csv


# e-features from EFel

sh e_features/dl_e_features.sh

