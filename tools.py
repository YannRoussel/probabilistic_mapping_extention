import umap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import scanpy.external as sce

from gtfparse import read_gtf

from scipy import sparse

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.ensemble import RandomForestClassifier

def plot_umap(X, labels, title):
    """

    :param X:
    :param labels:
    :param title:
    :return:
    """
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(X)

    for lbl in labels.unique():
        msk = labels == lbl
        plt.scatter(embedding[:, 0][msk], embedding[:, 1][msk], alpha=0.5)
    plt.legend(labels.unique())
    plt.xlabel("UMAP_1")
    plt.ylabel("UMAP_2")
    plt.title(title)

    return


def count_elements(array):
    """

    :param array:
    :return:
    """
    cts = {}
    for x in np.unique(array):
        msk_ = [y == x for y in array]
        cts[x] = len(array[msk_])

    return pd.DataFrame.from_dict(cts, orient='index', columns=["counts"])

def form_feature_based_clusters(features, n_clusters=20):
    """

    :param features:
    :param n_clusters:
    :return:
    """
    features_norm = StandardScaler().fit_transform(features)
    clustering = AgglomerativeClustering(n_clusters).fit(features_norm)
    clstr_labels = clustering.labels_

    return features_norm, clstr_labels

def create_map_(labels_true, labels):
    
    cts_list = []
    for true_lbl in np.unique(labels_true):
        msk_ = np.asarray([lbl == true_lbl for lbl in labels_true])

        cts_list.append(count_elements(labels[msk_]).rename({"counts":true_lbl}, axis=1))

    map_ = pd.concat(cts_list, axis=1)
    
    return map_

def select_best_features(data, labels_type, labels_region, threshold):
    
    mask_features = {}
    for title, labels in zip(["region", "type"], [labels_region, labels_type]):
        clf = RandomForestClassifier()
        clf.fit(data, labels)
        values = clf.feature_importances_
        mask_features[title] = [x > np.quantile(values, threshold) for x in values]
        
    selected_features = []
    for x in data.columns[mask_features["type"]]:
        if x not in data.columns[mask_features["region"]]:
            selected_features.append(x)
            
    features = {
        "type" : data.columns[mask_features["type"]], 
        "region" : data.columns[mask_features["region"]],
        "selected": np.asarray(selected_features)
    }
    
    return features


def create_count_maps(data, labels_type, labels_region, features, n_clusters):
    
    common_clusters = form_feature_based_clusters(data[features["selected"]], n_clusters=n_clusters)[1]
    
    p_map_counts = {}
    for region in np.unique(labels_region):
        msk_reg = np.asarray([x==region for x in labels_region])
        p_map_counts[region] = create_map_(np.asarray(labels_type)[msk_reg], common_clusters[msk_reg])
        
    
    return p_map_counts

def create_probabilistic_map(p_map_counts):
    
    region_list = list(p_map_counts.keys())
    clstr_idx = p_map_counts[region_list[0]].index.intersection(p_map_counts[region_list[1]].index)
    
    for region in region_list:
        p_map_counts[region] = p_map_counts[region].reindex(clstr_idx)
              
    p_regA_counts = p_map_counts[region_list[0]].div(np.sum(p_map_counts[region_list[0]], axis=1), axis=0).fillna(0)
    p_counts_regB = p_map_counts[region_list[1]].div(np.sum(p_map_counts[region_list[1]], axis=0), axis=1).fillna(0)
    
    p_regA_regB = p_regA_counts.T @ p_counts_regB
    
    return p_regA_regB

def reindex_on_common_index(df_A, df_B):
    
    common_idx = df_A.index.intersection(df_B.index)
    df_A_ = df_A.reindex(common_idx)
    df_B_ = df_B.reindex(common_idx)
    
    return df_A_, df_B_

def best_features_based_p_maps(data, lbls_type, lbls_reg, threshold, n_clusters,
                               reg1, reg2):
    """
    """
    # select best m-features that explain mtypes labels
    features = select_best_features(data, lbls_type, lbls_reg, threshold=threshold)
    # create probabilistic map between m-types from different regions
    count_maps_m = create_count_maps(data, lbls_type, lbls_reg, features, n_clusters=n_clusters)
    # put regions in order
    count_maps_m_1_2 = {reg1: count_maps_m[reg1], reg2: count_maps_m[reg2]} 
    p_map_1_2 = create_probabilistic_map(count_maps_m_1_2)

    count_maps_m_2_1 = {reg2: count_maps_m[reg2], reg1: count_maps_m[reg1]} 
    p_map_2_1 = create_probabilistic_map(count_maps_m_2_1)
    
    return p_map_1_2, p_map_2_1, features

def geneSelection_Kobak(data, threshold=0, atleast=10,
                  yoffset=.02, xoffset=5, decay=1.5, n=None,
                  plot=True, markers=None, genes=None, figsize=(6, 3.5),
                  markeroffsets=None, labelsize=10, alpha=1, verbose=1):
    
    import seaborn as sns
    
    if sparse.issparse(data):
        zeroRate = 1 - np.squeeze(np.array((data > threshold).mean(axis=0)))
        A = data.multiply(data > threshold)
        A.data = np.log2(A.data)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        meanExpr[detected] = np.squeeze(np.array(A[:, detected].mean(axis=0))) / (1 - zeroRate[detected])
    else:
        zeroRate = 1 - np.mean(data > threshold, axis=0)
        meanExpr = np.zeros_like(zeroRate) * np.nan
        detected = zeroRate < 1
        mask = data[:, detected] > threshold
        logs = np.zeros_like(data[:, detected]) * np.nan
        logs[mask] = np.log2(data[:, detected][mask])
        meanExpr[detected] = np.nanmean(logs, axis=0)

    lowDetection = np.array(np.sum(data > threshold, axis=0)).squeeze() < atleast
    zeroRate[lowDetection] = np.nan
    meanExpr[lowDetection] = np.nan

    if n is not None:
        up = 10
        low = 0
        for t in range(100):
            nonan = ~np.isnan(zeroRate)
            selected = np.zeros_like(zeroRate).astype(bool)
            selected[nonan] = zeroRate[nonan] > np.exp(-decay * (meanExpr[nonan] - xoffset)) + yoffset
            if np.sum(selected) == n:
                break
            elif np.sum(selected) < n:
                up = xoffset
                xoffset = (xoffset + low) / 2
            else:
                low = xoffset
                xoffset = (xoffset + up) / 2
        if verbose > 0:
            print('Chosen offset: {:.2f}'.format(xoffset))
    else:
        nonan = ~np.isnan(zeroRate)
        selected = np.zeros_like(zeroRate).astype(bool)
        selected[nonan] = zeroRate[nonan] > np.exp(-decay * (meanExpr[nonan] - xoffset)) + yoffset

    if plot:
        if figsize is not None:
            plt.figure(figsize=figsize)
        plt.ylim([0, 1])
        if threshold > 0:
            plt.xlim([np.log2(threshold), np.ceil(np.nanmax(meanExpr))])
        else:
            plt.xlim([0, np.ceil(np.nanmax(meanExpr))])
        x = np.arange(plt.xlim()[0], plt.xlim()[1] + .1, .1)
        y = np.exp(-decay * (x - xoffset)) + yoffset
        if decay == 1:
            plt.text(.4, 0.2, '{} genes selected\ny = exp(-x+{:.2f})+{:.2f}'.format(np.sum(selected), xoffset, yoffset),
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)
        else:
            plt.text(.4, 0.2,
                     '{} genes selected\ny = exp(-{:.1f}*(x-{:.2f}))+{:.2f}'.format(np.sum(selected), decay, xoffset,
                                                                                    yoffset),
                     color='k', fontsize=labelsize, transform=plt.gca().transAxes)

        plt.plot(x, y, color=sns.color_palette()[1], linewidth=2)
        xy = np.concatenate((np.concatenate((x[:, None], y[:, None]), axis=1), np.array([[plt.xlim()[1], 1]])))
        t = plt.matplotlib.patches.Polygon(xy, color=sns.color_palette()[1], alpha=.4)
        plt.gca().add_patch(t)

        plt.scatter(meanExpr, zeroRate, s=1, alpha=alpha, rasterized=True)
        if threshold == 0:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of zero expression')
        else:
            plt.xlabel('Mean log2 nonzero expression')
            plt.ylabel('Frequency of near-zero expression')
        plt.tight_layout()

        if markers is not None and genes is not None:
            if markeroffsets is None:
                markeroffsets = [(0, 0) for g in markers]
            for num, g in enumerate(markers):
                i = np.where(genes == g)[0]
                plt.scatter(meanExpr[i], zeroRate[i], s=10, color='k')
                dx, dy = markeroffsets[num]
                plt.text(meanExpr[i] + dx + .1, zeroRate[i] + dy, g, color='k', fontsize=labelsize)

    return selected