import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from tools import *
from sscx_cellwise_petilla_etypes import *


def load_labels_sscx_sc(petilla_labels, sscx_sc_e_features_df):
    """Load SSCX single-cell labels."""
    labels_sscx_sc = np.asarray([petilla_labels[x] for x in sscx_sc_e_features_df.index])
    print(f"Number of SSCX single-cell labels: {len(labels_sscx_sc)}")
    return labels_sscx_sc


def convert_json_to_pd_e_features(path):
    """Convert JSON files to a Pandas DataFrame of e-features."""
    e_features_collection = []

    for e_type in os.listdir(path):
        if e_type.endswith(".json"):
            with open(os.path.join(path, e_type), 'r') as json_file:
                data_dict = json.load(json_file)

            e_names, e_vals, e_stds = [], [], []

            for d in data_dict["efeatures"]:
                prot_name = "_".join(["Step", d["protocol_name"].split("_")[-1]]) if any(
                    x in d["protocol_name"] for x in ["IDthresh", "IDrest"]) else d["protocol_name"]
                e_name = "|".join([prot_name, d["efel_feature_name"]])

                if e_name not in e_names:
                    e_names.append(e_name)
                    e_vals.append(d["mean"])
                    e_stds.append(d["original_std"])

            e_features_collection.append(pd.DataFrame(e_vals, index=e_names, columns=[e_type.split(".")[0]]).T)

    e_features_df = pd.concat(e_features_collection, axis=0)
    return e_features_df


def plot_map(data, title, filename):
    """Plot and save heatmap of data."""
    plt.figure(figsize=(25,10))
    heatmap = plt.pcolor(data.values, cmap="jet")

    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            if 0 < data.values[y, x] < 1.:
                plt.text(x + 0.5, y + 0.5, int(100 * data.values[y, x]),
                         horizontalalignment='center',
                         verticalalignment='center')

    plt.colorbar(heatmap)
    plt.xticks(np.arange(len(data.columns)) + .5, data.columns, rotation=90)
    plt.yticks(np.arange(len(data.index)) + .5, data.index)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def inverse_dict(dict_):
    """Inverse a dictionary."""
    inv_dict_ = {e_type.replace("-", "_"): k for k, e_type_lst in dict_.items() for e_type in e_type_lst}
    return inv_dict_


def counts_elements(labels):
    """Count occurrences of each unique element in labels."""
    return {lbl: np.sum(labels == lbl) for lbl in np.unique(labels)}


def map_counts(label1, label2):
    """Map counts of elements in label2 for each unique element in label1."""
    cts_ = []
    for lbl1 in np.unique(label1):
        mask = (label1 == lbl1)
        counts = counts_elements(label2[mask])
        cts_.append(pd.DataFrame(counts.values(), index=counts.keys(), columns=[lbl1]))

    return pd.concat(cts_, axis=1)


def clean_dataframes(shared_columns, *dfs):
    """Clean dataframes to only include shared columns and drop NaN values."""
    cleaned_dfs = []
    for df in dfs:
        df_cleaned = df[shared_columns].dropna(axis=0, how="any")
        cleaned_dfs.append(df_cleaned)
        print(f"Cleaned dataframe shape: {np.shape(df_cleaned)}")
    return cleaned_dfs


def load_e_features(paths):
    """Load and clean e-features dataframes from given paths."""
    e_features_dfs = {name: convert_json_to_pd_e_features(path).dropna(axis=1, how="any")
                      for name, path in paths.items()}
    for name, df in e_features_dfs.items():
        if name == "gouwens":
            e_features_dfs[name] = e_features_dfs[name].rename(
                { x : int(x.split("_")[-1]) for x in e_features_dfs[name].index},
                axis=0
                )
        print(f"{name} dataframe shape: {np.shape(df)}")
    return e_features_dfs


def filter_labels_and_clean_dataframes(scala_34_e_features_df, scala_rt_e_features_df, sscx_e_features_df, sscx_sc_e_features_df, gouwens_e_features_df):
    """Filter labels and clean dataframes to only include shared columns and drop NaN values."""
    shared_columns = scala_34_e_features_df.columns.intersection(scala_rt_e_features_df.columns).intersection(
        sscx_e_features_df.columns).intersection(sscx_sc_e_features_df.columns).intersection(gouwens_e_features_df.columns)
    print(f"Number of shared columns: {len(shared_columns)}")

    cleaned_dfs = clean_dataframes(shared_columns, scala_34_e_features_df, scala_rt_e_features_df, sscx_e_features_df,
                                   sscx_sc_e_features_df, gouwens_e_features_df)
    return cleaned_dfs, shared_columns


def load_and_filter_metadata(metadata_paths, e_features_dfs):
    """Load metadata and filter based on quality."""

    metadata = {}
    for name, path in metadata_paths.items():
        if name == "gouwens":
            metadata[name] = pd.read_csv(path, index_col=["ephys_session_id"])
        else:
            metadata[name] = pd.read_csv(path, index_col=["Cell"])

    labels_scala_34 = metadata["scala_34"]['RNA family'].reindex(
        [str(x.split(".")[0].replace("cell_", "")) for x in e_features_dfs["scala_34"].index]
    )
    mask_clean_scala_34 = ["low quality" not in x for x in labels_scala_34]
    e_features_dfs["scala_34"] = e_features_dfs["scala_34"][mask_clean_scala_34]
    labels_scala_34 = labels_scala_34[mask_clean_scala_34]

    labels_scala_rt = metadata["scala_rt"]['RNA family'].reindex(
        [str(x.split(".")[0].replace("cell_", "")) for x in e_features_dfs["scala_rt"].index]
    )
    mask_clean_scala_rt = ["low quality" not in x for x in labels_scala_rt]
    e_features_dfs["scala_rt"] = e_features_dfs["scala_rt"][mask_clean_scala_rt]
    labels_scala_rt = labels_scala_rt[mask_clean_scala_rt]

    metadata_gouwens = pd.read_csv(metadata_paths["gouwens"], index_col="ephys_session_id")
    labels_gouwens = np.asarray([str(x) for x in metadata_gouwens["corresponding_AIT2.3.1_alias"][e_features_dfs["gouwens"].index]])
    mask_clean_gouwens = ["nan" not in x for x in labels_gouwens]
    e_features_dfs["gouwens"] = e_features_dfs["gouwens"][mask_clean_gouwens].dropna(axis=0, how="any")
    labels_gouwens = np.asarray([str(x).split(" ")[0] for x in metadata_gouwens["corresponding_AIT2.3.1_alias"][e_features_dfs["gouwens"].index]])

    return labels_scala_34, labels_scala_rt, labels_gouwens, e_features_dfs


def main():
    PATH_SCALA_34 = "../experimental_data/e_features/scala_34/"
    PATH_SCALA_RT = "../experimental_data/e_features/scala_roomtemp/"
    PATH_GOUWENS = "../experimental_data/e_features/allen_34_VCtx_partial/"
    PATH_SSCX = "../experimental_data/e_features/sscx/"
    PATH_SSCX_single_cell = "../experimental_data/e_features/sscx_cellwise/"

    paths = {
        "scala_34": PATH_SCALA_34,
        "scala_rt": PATH_SCALA_RT,
        "gouwens": PATH_GOUWENS,
        "sscx": PATH_SSCX,
        "sscx_single_cell": PATH_SSCX_single_cell
    }

    e_features_dfs = load_e_features(paths)
    cleaned_dfs, shared_columns = filter_labels_and_clean_dataframes(*e_features_dfs.values())
    e_features_cleaned_dfs = {x : y for (x, y) in zip(paths.keys(), cleaned_dfs)}

    metadata_paths = {
        "scala_34": "../experimental_data/e_features/m1_patchseq_phys_temp_meta_data.csv",
        "scala_rt": "../experimental_data/scala/m1_patchseq_meta_data.csv",
        "gouwens": "../experimental_data/gouwens/20200625_patchseq_metadata_mouse.csv"
    }

    labels_scala_34, labels_scala_rt, labels_gouwens, cleaned_dfs = load_and_filter_metadata(metadata_paths, e_features_cleaned_dfs)
    # print("labels", labels_scala_34, labels_scala_rt, labels_gouwens)
    sscx_sc_e_features_df, scala_rt_e_features_df = cleaned_dfs["sscx_single_cell"], cleaned_dfs["scala_rt"]

    petilla_labels = inverse_dict(cells_per_etype)
    labels_sscx_sc = load_labels_sscx_sc(petilla_labels, sscx_sc_e_features_df)

    scala_rt_yao_types = pd.read_csv("Scala_aligned_labels.csv", index_col=0)

    mask_exc = np.asarray(["Glut" in x for x in scala_rt_yao_types["AIBS lvl 4"]]) & np.asarray(
        ["excitatory" in y for y in scala_rt_yao_types["Scala superfamily"]])
    mask_inh = np.asarray(["Gaba" in x for x in scala_rt_yao_types["AIBS lvl 4"]]) & np.asarray(
        ["inhibitory" in y for y in scala_rt_yao_types["Scala superfamily"]])

    map_exc = map_counts(scala_rt_yao_types[mask_exc]["AIBS lvl 3"], scala_rt_yao_types[mask_exc]["Scala type"])
    map_inh = map_counts(scala_rt_yao_types[mask_inh]["AIBS lvl 3"], scala_rt_yao_types[mask_inh]["Scala type"])

    yao_types_kept = list(map_exc.columns[np.sum(map_exc, axis=0) >= 3]) + list(map_inh.columns[np.sum(map_inh, axis=0) >= 3])
    print(f"Yao types kept: {yao_types_kept}")

    mask_next_lvl_exc = np.asarray([x in yao_types_kept for x in scala_rt_yao_types["AIBS lvl 3"]]) & mask_exc
    mask_next_lvl_inh = np.asarray([x in yao_types_kept for x in scala_rt_yao_types["AIBS lvl 3"]]) & mask_inh

    scala_rt_yao_types_cleaned = pd.concat([scala_rt_yao_types[mask_next_lvl_exc], scala_rt_yao_types[mask_next_lvl_inh]], axis=0)

    scala_rt_e_features_df = scala_rt_e_features_df.rename({x: x.replace("cell_", "") for x in scala_rt_e_features_df.index})

    scaler_scala = StandardScaler()
    scaler_bbp = StandardScaler()

    e_index = scala_rt_yao_types_cleaned.index.intersection(scala_rt_e_features_df.index)

    scala_rt_e_features_df.reindex(e_index).to_csv("scala_rt_e_features_aligned_to_yao_t_types.csv")

    combined_e_feature_df = pd.concat([
        pd.DataFrame(scaler_scala.fit_transform(scala_rt_e_features_df.reindex(e_index)), index=e_index,
                     columns=scala_rt_e_features_df.columns),
        pd.DataFrame(scaler_bbp.fit_transform(sscx_sc_e_features_df), index=sscx_sc_e_features_df.index,
                     columns=sscx_sc_e_features_df.columns)
    ])

    scala_used_labels = scala_rt_yao_types_cleaned["AIBS lvl 0"][[x.replace("cell_", "") for x in e_index]]
    scala_used_labels = np.asarray(scala_used_labels)

    msk_scala_e_inh = np.asarray(["Gaba" in x for x in scala_used_labels])
    msk_scala_e_exc = np.asarray(["Glut" in x for x in scala_used_labels])

    combined_labels = list(scala_used_labels) + list(labels_sscx_sc)
    combined_labels = np.asarray(combined_labels)
    combined_dset_labels = ["scala_rt"] * len(scala_used_labels) + ["sscx_sc"] * len(labels_sscx_sc)

    p_e_t, e_features = best_features_based_p_maps(combined_e_feature_df, combined_labels, combined_dset_labels, 
                                                   threshold=.9, n_clusters=4,
                                                   reg1="sscx_sc", reg2="scala_rt")
    
    # Identify columns that contain 'Glut' in their names
    glut_columns = [col for col in p_e_t.columns if 'Glut' in col]

    # Replace values in 'Glut' columns
    for col in glut_columns:
        p_e_t[col] = p_e_t.index.to_series().apply(lambda x: 1. if x == 'cADpyr' else 0.)

    # Identify columns that contain 'Gaba' in their names
    gaba_columns = [col for col in p_e_t.columns if 'Gaba' in col]

    # Update values in 'Gaba' columns
    for col in gaba_columns:
        p_e_t.at['cAC', col] = p_e_t.at['cAC', col] + p_e_t.at['cADpyr', col]
        p_e_t.at['cADpyr', col] = 0

    p_e_t.to_csv("p_e_t.csv")
    plot_map(p_e_t, "P(sscx|rt) yao", "p_sscx_rt_yao.png")

    canonical_types = pd.read_csv("/Users/yroussel/Documents/data/Patch Seq data/test_cluster.csv", index_col=0, header=None)
    canonical_types = canonical_types.rename({i: x for i, x in enumerate(["name", "class", "dendrites", "axon"])}, axis=1)
    canonical_types = canonical_types.rename({x: x.split(".")[0] for x in canonical_types.index}, axis=0)

    gouwens_idx_dict = {x: x.split("_")[0] for x in canonical_types.index if "sample" not in x}
    canonical_types = canonical_types.rename(gouwens_idx_dict, axis=0)

    dict_mclass = {0.0: "IN", 1.0: "PC"}
    dict_dendclass = {x: "dend_" + str(int(x)) for x in np.unique(canonical_types["dendrites"])}
    dict_axclass = {x: "ax_" + str(int(x)) for x in np.unique(canonical_types["axon"])}
    dict_mlabels = {"class": dict_mclass, "dendrites": dict_dendclass, "axon": dict_axclass}

    m_index = scala_rt_yao_types_cleaned.index.intersection(canonical_types.index)
    mt_labels = pd.concat([scala_rt_yao_types_cleaned.T[m_index].T, canonical_types.T[m_index].T], axis=1)
    mt_labels.to_csv("t_and_m_types.csv")

    m_labels = ["_".join([dict_mlabels["class"][mt_labels["class"][x]],
                          dict_mlabels["dendrites"][mt_labels["dendrites"][x]],
                          dict_mlabels["axon"][mt_labels["axon"][x]]
                          ]) for x in mt_labels.index]
    t_labels = mt_labels["AIBS lvl 0"]

    m_labels = np.asarray(m_labels)
    t_labels = np.asarray(t_labels)

    cts_m_t = map_counts(m_labels, t_labels)
    p_t_m = cts_m_t.div(np.sum(cts_m_t, axis=0), axis=1)
    p_m_t = cts_m_t.div(np.sum(cts_m_t, axis=1), axis=0)

    plot_map(p_m_t, "P(canonical_m|yao_t)", "p_canonical_m_yao_t.png")
    p_m_t.to_csv("p_m_t.csv")

    p_e_t = p_e_t.T
    p_m_t = cts_m_t.div(np.sum(cts_m_t, axis=1), axis=0)

    index_union = p_e_t.index.union(p_m_t.index).unique()

    p_e_t = p_e_t.reindex(index_union)
    p_m_t = p_m_t.reindex(index_union)

    df_coll = []
    for col in p_e_t.columns:
        p_ = p_m_t.multiply(p_e_t[col], axis=0)
        df_coll.append(p_.rename({x: "|".join([x, col]) for x in p_m_t.columns}, axis=1))

    p_me_t = pd.concat(df_coll, axis=1)
    p_me_t = p_me_t.rename({x: x.replace(" ", "_") for x in p_me_t.index}, axis=0)

    plot_map(p_me_t, "P(bbp_e__canonical_m|yao_t)", "p_bbp_e_canonical_m_yao_t.png")
    p_me_t.to_csv("p_me_t.csv")


if __name__ == "__main__":
    main()
