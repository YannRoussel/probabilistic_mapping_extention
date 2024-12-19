#densities override
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt

# Load data
# Replace with actual data loading
t_types_expression = pd.read_csv('t_types_expression.csv', index_col=0)
t_types_densities = pd.read_csv('t_types_densities.csv', index_col=0)
me_types_densities = pd.read_csv('me_types_densities.csv', index_col=0)

# Step 1: Hierarchical Clustering of T-types
distance_matrix = pdist(t_types_expression.values, metric='euclidean')
linkage_matrix = linkage(distance_matrix, method='ward')

plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=t_types_expression.index)
plt.title('Hierarchical Clustering Dendrogram of T-types')
plt.xlabel('T-types')
plt.ylabel('Distance')
plt.show()

# Step 2: Optimization
def compute_similarity(grouped_t_types_densities, me_types_densities):
    grouped_t_types_normalized = grouped_t_types_densities.div(grouped_t_types_densities.sum(axis=1), axis=0)
    me_types_normalized = me_types_densities.div(me_types_densities.sum(axis=1), axis=0)
    distance_matrix = pairwise_distances(grouped_t_types_normalized, me_types_normalized, metric='euclidean')
    return distance_matrix

def get_grouped_densities(t_types_densities, clusters):
    grouped_densities = pd.DataFrame(index=t_types_densities.index)
    for cluster_id in np.unique(clusters):
        grouped_densities[f'group_{cluster_id}'] = t_types_densities.loc[:, clusters == cluster_id].sum(axis=1)
    return grouped_densities

def optimize_grouping(t_types_densities, me_types_densities, linkage_matrix, num_clusters):
    clusters = fcluster(linkage_matrix, num_clusters, criterion='maxclust')
    grouped_t_types_densities = get_grouped_densities(t_types_densities, clusters)
    similarity = compute_similarity(grouped_t_types_densities, me_types_densities)
    return clusters, similarity

num_clusters = 10  # Adjust based on desired level of grouping
clusters, similarity = optimize_grouping(t_types_densities, me_types_densities, linkage_matrix, num_clusters)

# Step 3: Mapping
def map_t_types_to_me_types(similarity_matrix):
    mapping = {}
    for i in range(similarity_matrix.shape[0]):
        closest_me_type = np.argmin(similarity_matrix[i])
        mapping[f'group_{i+1}'] = me_types_densities.columns[closest_me_type]
    return mapping

mapping = map_t_types_to_me_types(similarity)
print(mapping)
