import mygene
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import sparse
import umap
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import scipy.spatial as spatial
from scipy.spatial.distance import pdist, squareform

def count_elements(array):
    array = np.asarray(array)
    cts = {}
    for x in np.unique(array):
        msk_ = [y == x for y in array]
        cts[x] = len(array[msk_])
    return cts

def geneSelection(data, threshold=0, atleast=10, 
                  yoffset=.02, xoffset=5, decay=1.5, n=None, 
                  plot=True, markers=None, genes=None, figsize=(6,3.5),
                  markeroffsets=None, labelsize=10, alpha=1, verbose=1):
    
    if sparse.issparse(data):
        zeroRate = 1 - np.squeeze(np.array((data>threshold).mean(axis=0)))
        A = data.multiply(data>threshold)
        A.data = np.log2(A.data)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        meanExpr[detected] = np.squeeze(np.array(A[:,detected].mean(axis=0))) / (1-zeroRate[detected])
    else:
        zeroRate = 1 - np.mean(data>threshold, axis=0)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        mask = data[:,detected]>threshold
        logs = np.zeros_like(data[:,detected]) * np.nan
        logs[mask] = np.log2(data[:,detected][mask])
        meanExpr[detected] = np.nanmean(logs, axis=0)

    lowDetection = np.array(np.sum(data>threshold, axis=0)).squeeze() < atleast
    zeroRate[lowDetection] = np.nan
    meanExpr[lowDetection] = np.nan
            
    if n is not None:
        up = 10
        low = 0
        for t in range(100):
            nonan = ~np.isnan(zeroRate)
            selected = np.zeros_like(zeroRate).astype(bool)
            selected[nonan] = zeroRate[nonan] > np.exp(-decay*(meanExpr[nonan] - xoffset)) + yoffset
            if np.sum(selected) == n:
                break
            elif np.sum(selected) < n:
                up = xoffset
                xoffset = (xoffset + low)/2
            else:
                low = xoffset
                xoffset = (xoffset + up)/2
        if verbose>0:
            print('Chosen offset: {:.2f}'.format(xoffset))
    else:
        nonan = ~np.isnan(zeroRate)
        selected = np.zeros_like(zeroRate).astype(bool)
        selected[nonan] = zeroRate[nonan] > np.exp(-decay*(meanExpr[nonan] - xoffset)) + yoffset
                
    if plot:
        if figsize is not None:
            plt.figure(figsize=figsize)
        plt.ylim([0, 1])
        if threshold>0:
            plt.xlim([np.log2(threshold), np.ceil(np.nanmax(meanExpr))])
        else:
            plt.xlim([0, np.ceil(np.nanmax(meanExpr))])
        x = np.arange(plt.xlim()[0], plt.xlim()[1]+.1,.1)
        y = np.exp(-decay*(x - xoffset)) + yoffset
        if decay==1:
            plt.text(.4, 0.2, '{} genes selected\ny = exp(-x+{:.2f})+{:.2f}'.format(np.sum(selected),xoffset, yoffset), 
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)
        else:
            plt.text(.4, 0.2, '{} genes selected\ny = exp(-{:.1f}*(x-{:.2f}))+{:.2f}'.format(np.sum(selected),decay,xoffset, yoffset), 
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)

        plt.plot(x, y, color=sns.color_palette()[1], linewidth=2)
        xy = np.concatenate((np.concatenate((x[:,None],y[:,None]),axis=1), np.array([[plt.xlim()[1], 1]])))
        t = plt.matplotlib.patches.Polygon(xy, color=sns.color_palette()[1], alpha=.4)
        plt.gca().add_patch(t)
        
        plt.scatter(meanExpr, zeroRate, s=1, alpha=alpha, rasterized=True)
        if threshold==0:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of zero expression')
        else:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of near-zero expression')
        plt.tight_layout()
        
        if markers is not None and genes is not None:
            if markeroffsets is None:
                markeroffsets = [(0, 0) for g in markers]
            for num,g in enumerate(markers):
                i = np.where(genes==g)[0]
                plt.scatter(meanExpr[i], zeroRate[i], s=10, color='k')
                dx, dy = markeroffsets[num]
                plt.text(meanExpr[i]+dx+.1, zeroRate[i]+dy, g, color='k', fontsize=labelsize)
    
    return selected

def preprocess_RNA(RNA_, kept_genes , volume_correction=False, volume_correction_df=None):
    X = (RNA_[kept_genes].div(np.sum(RNA_, axis=1), axis=0) * 1e+6)
    X = np.array(X)
    X = np.log10(X + 1)
    if volume_correction:
        X = X / volume_correction_df.reindex(RNA_.index).values
    X = X - np.mean(X, axis=0)
    X = X / np.std(X, axis=0)
    X = pd.DataFrame(X, index=RNA_.index, columns=kept_genes)
    
    return X
    
def plot_umap(X, labels, title):

    reducer = umap.UMAP()
    embedding = reducer.fit_transform(X)

    for lbl in np.unique(labels):
        msk = labels == lbl
        plt.plot(embedding[:,0][msk], embedding[:,1][msk], "o", alpha=0.5, ms=8)
    plt.legend(np.unique(labels))
    plt.xlabel("UMAP_1")
    plt.ylabel("UMAP_2")
    plt.title(title)
#     plt.show()
    return

def compute_distance_matrix(X):
    """
    Compute distance matrix between each pair of samples in matrix X.

    Parameters:
    X : array-like, shape (n_samples, n_features)
        Input matrix of samples.

    Returns:
    distance_matrix : array, shape (n_samples, n_samples)
        Distance matrix between each pair of samples.
    """
    # Compute pairwise distances between samples
    pairwise_distances = pdist(X)
    
    # Convert pairwise distances to a square distance matrix
    distance_matrix = squareform(pairwise_distances)
    
    return distance_matrix

def find_closest_row(X, Y):
  """
  This function takes as input two matrices X, Y of shapes (nx, nf) and (ny, nf), respectively, and returns the row index of X that is the closest to each row of Y.

  Args:
    X: A matrix of shape (nx, nf).
    Y: A matrix of shape (ny, nf).

  Returns:
    A list of row indices of X that are the closest to each row of Y.
  """

  closest_row_indices = []
  min_row_distances = []
  ref_dist = np.median(pdist(X))

  for i in range(len(Y)):
    row_distances = []
    for j in range(len(X)):
      # row_distances.append(np.linalg.norm(X[j] - Y[i]))
      row_distances.append(spatial.distance.euclidean(X[j], Y[i]))
    
    # if min(row_distances) < ref_dist:
    closest_row_indices.append(row_distances.index(min(row_distances)))
    # else:
    #   print("no match")
    #   closest_row_indices.append(-1)

  return closest_row_indices
       
def map_counts(label1, label2):
    cts_ = []
    for lbl1 in np.unique(label1):
        msk_ = [l==lbl1 for l in label1]
        cts = count_elements(label2[msk_])
        cts_.append(pd.DataFrame(cts.values(), index=cts.keys(), columns=[lbl1]))

    return pd.concat(cts_, axis=1)

def plot_map(data, title):
    
    heatmap = plt.pcolor(data.values, cmap="jet")

    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            if 0<data.values[y, x]<1.:
                plt.text(x + 0.5, y + 0.5, int(100*data.values[y, x]),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontsize=8
                         )

    plt.colorbar(heatmap)
    plt.xticks(np.arange(len(data.columns)) + .5, data.columns, rotation=90)
    plt.yticks(np.arange(len(data.index)) + .5, data.index)
    plt.title(title)
#     plt.show()
    
    return

def convert_gene_tags_to_labels(dataframe):
    """
    Convert gene tags to gene labels using MyGeneInfo API.

    :param dataframe: Pandas DataFrame containing gene tags as column names.
    :return: Pandas DataFrame with gene labels replacing gene tags.
    """
    mg = mygene.MyGeneInfo()
    ens = dataframe.columns.to_list()
    ginfo = mg.querymany(ens, scopes='ensembl.gene')

    dict_id_to_symbol = {}

    for g in ginfo:
        try:
            dict_id_to_symbol[g["query"]] = g["symbol"]
        except:
            dict_id_to_symbol[g["query"]] = g["query"]

    return dataframe.rename(dict_id_to_symbol, axis=1)

def merge_and_label_RNA_data(ref_data, sample_data, ref_label="MERFISH", sample_label="Scala"):
    """
    Merge RNA data from reference and sample datasets while labeling each sample.

    :param ref_data: Pandas DataFrame containing RNA data from the reference dataset.
    :param sample_data: Pandas DataFrame containing RNA data from the sample dataset.
    :param ref_label: Label for the reference dataset (default is "MERFISH").
    :param sample_label: Label for the sample dataset (default is "Scala").
    :return: Merged Pandas DataFrame with labeled samples.
    """
    # Only take common genes between the reference and sample data
    common_genes = ref_data.columns.intersection(sample_data.columns)
    ref_data = ref_data.reindex(common_genes, axis=1)
    sample_data = sample_data.reindex(common_genes, axis=1)

    # Merge RNA data
    merged_data = pd.concat([ref_data, sample_data], axis=0)

    # Create labels for the samples
    labels = [ref_label] * len(ref_data) + [sample_label] * len(sample_data)
    labels = np.asarray(labels)

    return merged_data, labels

def merge_datasets(yao_df, Scala_meta_data, RNA, unique_genes):
    """
    Merge MERFISH and Scala datasets.

    :param yao_df: Pandas DataFrame containing MERFISH data.
    :param Scala_meta_data: Pandas DataFrame containing Scala meta data.
    :param RNA: Pandas DataFrame containing Scala RNA data.
    :param unique_genes: List of unique genes.
    :return: AnnData object containing merged datasets.
    """
    yao_df = yao_df[unique_genes].dropna(how="any", axis=1)

    common_genes = yao_df.columns.intersection(RNA.columns)

    scala_obs = pd.DataFrame(Scala_meta_data.reindex(RNA.index)["RNA family"].to_list(), columns=["dataset"])

    adata_scala = ad.AnnData(RNA.reindex(common_genes, axis=1).values, obs=scala_obs, 
                             var=pd.DataFrame(common_genes, columns=["gene_symbol"]), dtype='float64')

    merfish_obs = pd.DataFrame(["MERFISH"] * len(yao_df.index), columns=["dataset"])

    adata_merfish = ad.AnnData(yao_df.reindex(common_genes, axis=1).values, obs=merfish_obs, 
                             var=pd.DataFrame(common_genes, columns=["gene_symbol"]), dtype='float64')

    var_names = adata_merfish.var_names.intersection(adata_scala.var_names)
    adata_merfish = adata_merfish[:, var_names]
    adata_scala = adata_scala[:, var_names]

    sc.pp.pca(adata_merfish)
    sc.pp.neighbors(adata_merfish)
    sc.tl.umap(adata_merfish)

    sc.tl.ingest(adata_scala, adata_merfish, obs='dataset')

    adata_concat = adata_merfish.concatenate(adata_scala, batch_categories=['merfish', 'scala'])
    sc.pl.umap(adata_concat, color='batch')

    return adata_scala, adata_merfish, adata_concat

def compute_aligned_labels(adata_merfish, adata_scala, yao_df, RNA, msk_ctx, Scala_meta_data):
    """
    Compute aligned labels between MERFISH and Scala datasets.

    :param adata_merfish: AnnData object containing MERFISH data.
    :param adata_scala: AnnData object containing Scala data.
    :param yao_df: DataFrame containing preprocessed MERFISH data.
    :param RNA: DataFrame containing RNA data.
    :param msk_ctx: Mask indicating the context for alignment.
    :param Scala_meta_data: DataFrame containing Scala meta data.
    :return: DataFrame containing aligned labels.
    """
    indices = find_closest_row(adata_merfish.obsm["X_pca"][msk_ctx], adata_scala.obsm["X_pca"])

    scala_aligned_names = []
    s_names = []
    for n, idc in enumerate(indices):
        idc_ = yao_df[msk_ctx].index[idc]
        lbl_lvl0 = idc_
        lbl_other = "_".join(lbl_lvl0.split("_")[0].split(" ")[1:])
        lbl_lvl4 = lbl_other.split("_")[-1]
        lbl_lvl3 = "_".join(lbl_other.split("_")[-2:])
        lbl_lvl2 = "_".join(lbl_other.split("_")[-3:])
        lbl_lvl1 = "_".join(lbl_other.split("_")[-4:])

        table = [lbl_lvl0, lbl_lvl1, lbl_lvl2, lbl_lvl3, lbl_lvl4] + [Scala_meta_data.reindex(RNA.index)["RNA family"].to_list()[n]] + [Scala_meta_data.reindex(RNA.index)["RNA type"].to_list()[n]]

        scala_aligned_names.append(np.asarray(table))
        s_names.append(Scala_meta_data.reindex(RNA.index)["RNA family"].index[n]) 

    alignement = pd.DataFrame(
    np.asarray(scala_aligned_names), 
    columns = ["AIBS lvl 0", "AIBS lvl 1", "AIBS lvl 2", "AIBS lvl 3", "AIBS lvl 4", "Scala family", "Scala type"],
    index = s_names
    )

    dict_scala_ei = {
    "Lamp5": "inhibitory", "Vip": "inhibitory",  "Sncg": "inhibitory", "Pvalb": "inhibitory", "Sst": "inhibitory",
    "CT": "excitatory", "ET": "excitatory", "IT": "excitatory", "NP": "excitatory",
                }
    column_to_add = pd.DataFrame([dict_scala_ei[x] for x in alignement["Scala family"]], index = alignement.index, columns=["Scala superfamily"])

    alignement_complete = pd.concat([alignement, column_to_add], axis=1)

    return alignement_complete

def compute_and_plot_maps(alignement_complete,mask, mask_name):
    map_ = map_counts(alignement_complete[mask]["AIBS lvl 3"], alignement_complete[mask]["Scala type"])
    yao_types_kept = map_.columns[np.sum(map_, axis=0) >= 3]
    mask_next_lvl = np.asarray([x in yao_types_kept for x in alignement_complete["AIBS lvl 3"]]) & mask
    map_next_lvl = map_counts(alignement_complete[mask_next_lvl]["AIBS lvl 1"], alignement_complete[mask_next_lvl]["Scala type"])

    plt.figure(figsize=(10, 10))
    plt.subplot(121)
    plot_map(map_.div(np.sum(map_, axis=1), axis=0), f"P(old|new) ({mask_name})")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")

    plt.subplot(122)
    plot_map(map_.div(np.sum(map_, axis=0), axis=1), f"P(new|old) ({mask_name})")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")

def plot_combined_maps(map_exc, map_inh):

    tick_label_fontsize = 8
    plt.figure(figsize=(11, 15))

    # Plot maps for mask_exc
    plt.subplot(221)
    plot_map(map_exc.div(np.sum(map_exc, axis=1), axis=0), "P(new|old) (Excitatory)")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")
    plt.xticks(fontsize=tick_label_fontsize)  # Adjust x-axis tick label fontsize
    plt.yticks(fontsize=tick_label_fontsize)

    plt.subplot(222)
    plot_map(map_exc.div(np.sum(map_exc, axis=0), axis=1), "P(old|new) (Excitatory)")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")
    plt.xticks(fontsize=tick_label_fontsize)  # Adjust x-axis tick label fontsize
    plt.yticks(fontsize=tick_label_fontsize)

    # Plot maps for mask_inh
    plt.subplot(223)
    plot_map(map_inh.div(np.sum(map_inh, axis=1), axis=0), "P(new|old) (Inhibitory)")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")
    plt.xticks(fontsize=tick_label_fontsize)  # Adjust x-axis tick label fontsize
    plt.yticks(fontsize=tick_label_fontsize)

    plt.subplot(224)
    plot_map(map_inh.div(np.sum(map_inh, axis=0), axis=1), "P(old|new) (Inhibitory)")
    plt.xlabel("New t-type")
    plt.ylabel("Old t-type")
    plt.xticks(fontsize=tick_label_fontsize)  # Adjust x-axis tick label fontsize
    plt.yticks(fontsize=tick_label_fontsize)

    plt.tight_layout()
    plt.savefig("maps_combined.png")  # Save the figure
    plt.show()

def main():

    # get reference t-types
    yao = pd.read_csv("../experimental_data/yao/all_clusters_trim_mean25_g32285ct5322.csv", index_col="gene_identifier").T
    yao = convert_gene_tags_to_labels(yao)
    # clean genes of ref t-types
    counts_duplicates = count_elements(yao.columns)

    unique_genes = []
    for gene in counts_duplicates.keys():
        if counts_duplicates[gene] < 2:
            unique_genes.append(gene)

    # get data to align to the reference data
    PATH = "../experimental_data/scala"

    Scala_meta_data = pd.read_csv(PATH + '/m1_patchseq_meta_data.csv', index_col='Cell')
    data_filter = np.asarray([x != "low quality" for x in Scala_meta_data["RNA family"]])
    Scala_meta_data = Scala_meta_data[data_filter]

    RNA_e = pd.read_csv(PATH + "/rna/m1_patchseq_exon_counts.csv", index_col=0).T
    RNA_i = pd.read_csv(PATH + "/rna/m1_patchseq_intron_counts.csv", index_col=0).T

    RNA = (RNA_e + RNA_i)
    RNA = RNA.reindex(Scala_meta_data.index)

    # add gene selection on RNA here (and merge with Gouwens)
    selected_genes_scala = geneSelection(RNA.values, yoffset=.04, xoffset=6.2, decay=1.2)
    # selected_genes_scala = geneSelection(RNA.values, yoffset=.04, xoffset=6, decay=1.2)
    selected_genes_scala = RNA.columns[selected_genes_scala]
    RNA = RNA[selected_genes_scala]

    adata_scala, adata_merfish, adata_concat = merge_datasets(yao, Scala_meta_data, RNA, unique_genes)

    msk_ctx = [("NN_" not in x) & ("IMN_" not in x) for x in yao.index]
    # msk_ctx = [("NoMask_" not in x) for x in yao.index]
    # msk_ctx = [("OB-STR-CTX" in x) | ("Pvalb" in x) for x in vals_preprocess_df.index]
    # msk_ctx = [(("CTX" in x) | ("Pvalb" in x) | ("Sst" in x) | ("Vip" in x) | ("Lamp5" in x) | ("Sncg" in x)) & ("OB-STR-CTX" not in x)
    #            for x in vals_preprocess_df.index]

    # assign the closest t-type of the reference dataset for each cells of the data to align

    alignement_complete = compute_aligned_labels(adata_merfish, adata_scala, yao, RNA, msk_ctx, Scala_meta_data)
    alignement_complete.to_csv("./Scala_aligned_labels.csv")

    # make some figures to illustrate

    mask_exc = np.asarray(["Glut" in x for x in alignement_complete["AIBS lvl 4"]]) & np.asarray(["excitatory" in y for y in alignement_complete["Scala superfamily"]])
    mask_inh = np.asarray(["Gaba" in x for x in alignement_complete["AIBS lvl 4"]]) & np.asarray(["inhibitory" in y for y in alignement_complete["Scala superfamily"]])

    map_exc = map_counts(alignement_complete[mask_exc]["AIBS lvl 3"], alignement_complete[mask_exc]["Scala type"])
    map_inh = map_counts(alignement_complete[mask_inh]["AIBS lvl 3"], alignement_complete[mask_inh]["Scala type"])

    compute_and_plot_maps(alignement_complete,mask_exc, "Excitatory")
    compute_and_plot_maps(alignement_complete, mask_inh, "Inhibitory")

    plot_combined_maps(map_exc, map_inh)
    plt.show()

if __name__ == "__main__":
    main()