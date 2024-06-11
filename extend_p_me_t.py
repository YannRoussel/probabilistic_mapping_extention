import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import MultiTaskLassoCV, LogisticRegressionCV
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, KFold
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

# Load the data
yao_RNA = pd.read_csv("/Users/yroussel/Documents/projects/t-type_alignment/all_clusters_trim_mean25_g32285ct5322.csv", index_col="gene_identifier").T
t_tpes_dict = {x: x.replace(" ", "_") for x in yao_RNA.index}
yao_RNA = yao_RNA.rename(t_tpes_dict, axis=0)

existing_p_me_t = pd.read_csv("./p_me_t.csv", index_col=0)
covered_t_types = existing_p_me_t.index

uncovered_t_types = [x for x in yao_RNA.index if x not in covered_t_types]
gene_expression_df = yao_RNA.T[covered_t_types].T
uncovered_t_types_df = yao_RNA.T[uncovered_t_types].T

# Gene selection functions
def select_genes_for_efeatures(gene_expression, e_features, method='lasso', alpha=0.01, max_iter=5000):
    if method == 'lasso':
        model = MultiTaskLassoCV(alphas=[alpha], cv=5, max_iter=max_iter).fit(gene_expression.values, e_features.values)
    elif method == 'random_forest':
        model = RandomForestRegressor(n_estimators=100, random_state=42).fit(gene_expression.values, e_features.values)
    else:
        raise ValueError("Method not supported")
    importance_scores = np.abs(model.coef_).mean(axis=0) if method == 'lasso' else model.feature_importances_
    return pd.DataFrame({'gene': gene_expression.columns, 'importance': importance_scores}).sort_values(by='importance', ascending=False)

def select_genes_for_mtypes(gene_expression, m_types, method='logistic_regression', cv=5, random_state=42, max_iter=100):
    if method == 'logistic_regression':
        if len(np.unique(m_types)) < 2:
            raise ValueError("Logistic regression requires at least two classes in the target variable.")
        model = LogisticRegressionCV(cv=cv, random_state=random_state, max_iter=max_iter).fit(gene_expression.values, m_types)
        importance_scores = np.abs(model.coef_).mean(axis=0)
    elif method == 'random_forest':
        model = RandomForestClassifier(n_estimators=100, random_state=random_state).fit(gene_expression.values, m_types)
        importance_scores = model.feature_importances_
    elif method == 'mutual_info':
        importance_scores = mutual_info_classif(gene_expression.values, m_types, random_state=random_state)
    else:
        raise ValueError("Method not supported")
    return pd.DataFrame({'gene': gene_expression.columns, 'importance': importance_scores}).sort_values(by='importance', ascending=False)

def combine_feature_selection(gene_importance_df_e, gene_importance_df_m, top_n=100):
    combined_df = pd.merge(gene_importance_df_e, gene_importance_df_m, on='gene', suffixes=('_e', '_m'))
    combined_df['normalized_importance_e'] = combined_df['importance_e'] / combined_df['importance_e'].max()
    combined_df['normalized_importance_m'] = combined_df['importance_m'] / combined_df['importance_m'].max()
    combined_df['combined_importance'] = combined_df['normalized_importance_e'] + combined_df['normalized_importance_m']
    return combined_df.sort_values(by='combined_importance', ascending=False).head(top_n)

# Load and preprocess additional data
scala_rt_e_features_df = pd.read_csv("scala_rt_e_features_aligned_to_yao_t_types.csv", index_col=0)
mt_labels = pd.read_csv("t_and_m_types.csv", index_col=0)

common_samples = scala_rt_e_features_df.index.intersection(mt_labels.index)
scala_rt_e_features_df, mt_labels = scala_rt_e_features_df.reindex(common_samples), mt_labels.reindex(common_samples)

dict_mclass = {0.0: "IN", 1.0: "PC"}
dict_dendclass = {x: "dend_" + str(int(x)) for x in np.unique(mt_labels["dendrites"])}
dict_axclass = {x: "ax_" + str(int(x)) for x in np.unique(mt_labels["axon"])}
dict_mlabels = {"class": dict_mclass, "dendrites": dict_dendclass, "axon": dict_axclass}

m_labels = ["_".join([dict_mlabels["class"][mt_labels["class"][x]], dict_mlabels["dendrites"][mt_labels["dendrites"][x]], dict_mlabels["axon"][mt_labels["axon"][x]]]) for x in mt_labels.index]
t_labels = mt_labels["AIBS lvl 0"]

gene_expression = pd.concat([yao_RNA.T[t_tpes_dict[t_lbl]].T for t_lbl in t_labels], axis=1).T
m_labels = np.asarray(m_labels)
t_labels = np.asarray(t_labels)

# Error metrics
def mean_squared_error(actual, predicted):
    return np.mean((actual - predicted) ** 2)

def mean_absolute_error(actual, predicted):
    return np.mean(np.abs(actual - predicted))

def kullback_leibler_divergence(p, q):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

# Evaluation
def evaluate_n_choices(covered_p_me_t, gene_expression_scaled_df, n_choices, m_labels=None):
    errors = {n: [] for n in n_choices}
    
    if m_labels is not None and len(np.unique(m_labels)) > 1:
        strat_kf = StratifiedKFold(n_splits=5)
        split = strat_kf.split(gene_expression_scaled_df, m_labels)
    else:
        kf = KFold(n_splits=5)
        split = kf.split(gene_expression_scaled_df)
    
    for train_index, test_index in split:
        train_data = gene_expression_scaled_df.iloc[train_index]
        test_data = gene_expression_scaled_df.iloc[test_index]
        actual_p_me_t = covered_p_me_t.iloc[test_index]
        
        train_distance_matrix = cdist(test_data, train_data, metric='euclidean')
        
        for n in n_choices:
            predicted_p_me_t = []
            for i, test_idx in enumerate(test_index):
                closest_indices = np.argsort(train_distance_matrix[i])[:n]
                closest_t_types = train_data.index[closest_indices]
                avg_probabilities = covered_p_me_t.loc[closest_t_types].mean(axis=0)
                predicted_p_me_t.append(avg_probabilities)
                
            predicted_p_me_t = pd.DataFrame(predicted_p_me_t, index=actual_p_me_t.index, columns=actual_p_me_t.columns)
            
            mse = mean_squared_error(actual_p_me_t, predicted_p_me_t)
            mae = mean_absolute_error(actual_p_me_t, predicted_p_me_t)
            kld = kullback_leibler_divergence(actual_p_me_t.values.flatten(), predicted_p_me_t.values.flatten())
            
            errors[n].append((mse, mae, kld))
    
    return errors

def analyze_errors(errors):
    avg_errors = {n: np.mean(errors[n], axis=0) for n in errors}
    best_n = min(avg_errors, key=lambda n: avg_errors[n][0])  # Using MSE as the primary criterion
    
    for n, error_values in avg_errors.items():
        print(f"N={n}: MSE={error_values[0]:.4f}, MAE={error_values[1]:.4f}, KLD={error_values[2]:.4f}")
    
    return best_n

# Benchmarking and plotting
def benchmark_and_plot(covered_p_me_t, gene_expression, e_features, m_labels, n_choices, top_genes_nums):
    all_errors = {}
    e_methods = ['lasso', 'random_forest']
    # e_methods = ['random_forest']
    m_methods = ['logistic_regression', 'random_forest', 'mutual_info']
    # m_methods = ['random_forest']
    
    for e_method in e_methods:
        for m_method in m_methods:
            try:
                gene_importance_e = select_genes_for_efeatures(gene_expression, e_features, method=e_method)
                gene_importance_m = select_genes_for_mtypes(gene_expression, m_labels, method=m_method)
            except ValueError as e:
                print(f"Skipping combination {e_method} and {m_method} due to error: {e}")
                continue

            top_genes = combine_feature_selection(gene_importance_e, gene_importance_m, top_n=max(top_genes_nums))
            top_gene_names = top_genes['gene'].values
            gene_expression_filtered = gene_expression[top_gene_names]

            # Scaling
            combined_df = pd.concat([gene_expression_df[top_gene_names], uncovered_t_types_df[top_gene_names]])
            scaler = StandardScaler()
            combined_scaled = scaler.fit_transform(combined_df)
            gene_expression_scaled_df = pd.DataFrame(combined_scaled[:len(gene_expression_df)], index=gene_expression_df.index, columns=top_gene_names)
            uncovered_t_types_scaled_df = pd.DataFrame(combined_scaled[len(gene_expression_df):], index=uncovered_t_types_df.index, columns=top_gene_names)

            for top_n in top_genes_nums:
                gene_expression_top_n = gene_expression_scaled_df.iloc[:, :top_n]
                # gene_expression_top_n = gene_expression_filtered.iloc[:, :top_n]
                errors = evaluate_n_choices(existing_p_me_t, gene_expression_top_n, n_choices, m_labels=None)
                all_errors[(e_method, m_method, top_n)] = errors
            
            fig, axs = plt.subplots(3, 1, figsize=(10, 15))
            metrics = ['MSE', 'MAE', 'KLD']
    
            for i, metric in enumerate(metrics):
                for top_n in top_genes_nums:
                    avg_errors = {n: np.mean(all_errors[(e_method, m_method, top_n)][n], axis=0) for n in n_choices}
                    axs[i].plot(n_choices, [avg_errors[n][i] for n in n_choices], label=f'Top {top_n} genes')
                
                axs[i].set_title(f'{metric} vs Number of Neighbors (N) for {e_method} and {m_method}')
                axs[i].set_xlabel('Number of Neighbors (N)')
                axs[i].set_ylabel(metric)
                axs[i].legend()
    
            plt.tight_layout()
            plt.savefig(f'benchmark_{e_method}_{m_method}.png')
            plt.show()

# Run benchmark
n_choices = [1, 2, 3, 5, 8, 10, 20]
top_genes_nums = [50, 100, 500, 1000, 10000]
benchmark_and_plot(existing_p_me_t, gene_expression, scala_rt_e_features_df, m_labels, n_choices, top_genes_nums)

# Extend probabilistic map using best approach
def extend_probabilistic_map(existing_p_me_t, uncovered_t_types_scaled_df, gene_expression_scaled_df, best_n):
    distance_matrix = cdist(uncovered_t_types_scaled_df, gene_expression_scaled_df, metric='euclidean')
    uncovered_p_me_t = pd.DataFrame(index=uncovered_t_types_scaled_df.index, columns=existing_p_me_t.columns)
    
    for i, uncovered_t_type in enumerate(uncovered_t_types_scaled_df.index):
        closest_indices = np.argsort(distance_matrix[i])[:best_n]
        closest_t_types = gene_expression_scaled_df.index[closest_indices]
        avg_probabilities = existing_p_me_t.loc[closest_t_types].mean(axis=0)
        uncovered_p_me_t.loc[uncovered_t_type] = avg_probabilities
    
    uncovered_p_me_t.fillna(0, inplace=True)
    updated_p_me_t = pd.concat([existing_p_me_t, uncovered_p_me_t], axis=0)
    return updated_p_me_t

gene_importance_e = select_genes_for_efeatures(gene_expression, scala_rt_e_features_df, method='lasso')
gene_importance_m = select_genes_for_mtypes(gene_expression, m_labels, method='random_forest')
top_genes = combine_feature_selection(gene_importance_e, gene_importance_m, top_n=100)
top_gene_names = top_genes['gene'].values
gene_expression_filtered = gene_expression[top_gene_names]

# Scaling
combined_df = pd.concat([gene_expression_df[top_gene_names], uncovered_t_types_df[top_gene_names]])
scaler = StandardScaler()
combined_scaled = scaler.fit_transform(combined_df)
gene_expression_scaled_df = pd.DataFrame(combined_scaled[:len(gene_expression_df)], index=gene_expression_df.index, columns=top_gene_names)
uncovered_t_types_scaled_df = pd.DataFrame(combined_scaled[len(gene_expression_df):], index=uncovered_t_types_df.index, columns=top_gene_names)


best_n = 8#analyze_errors(evaluate_n_choices(existing_p_me_t, gene_expression_scaled_df, n_choices, m_labels=None))
extended_p_me_t = extend_probabilistic_map(existing_p_me_t, uncovered_t_types_scaled_df, gene_expression_scaled_df, best_n)
extended_p_me_t.to_csv("extended_p_me_t.csv")
