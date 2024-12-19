import pandas as pd
import numpy as np

extended_p_me_t = pd.read_csv("../probabilistic_mapping_extention/extended_p_me_t.csv", index_col=0)

common_t_types = extended_p_me_t.index
msk_t_types = np.asarray([("IMN" not in x) & ("NN" not in x) for x in common_t_types])
common_t_types = common_t_types[msk_t_types]

extended_p_me_t = extended_p_me_t.reindex(common_t_types, axis=0)
extended_p_me_t = extended_p_me_t.div(np.sum(extended_p_me_t, axis=1), axis=0)

df_col = []
for t in extended_p_me_t.index:
    df_ = extended_p_me_t.T[t].T
    df_ = df_.rename({x : "|".join([x, t]) for x in extended_p_me_t.columns}, axis=0)
    df_col.append(df_)

p_met_t = pd.concat(df_col, axis=1)

p_met_t