#!/usr/bin/env python
# coding: utf-8

import os
import re
import pandas as pd
import matplotlib.pyplot as plt

class NeuronLayerAnalyzer:
    def __init__(self, csv_path="data.csv", output_dir="results"):
        self.csv_path = csv_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        # Load CSV
        self.m_data_df = pd.read_csv(self.csv_path, index_col=0)

        # Prepare both collapsed versions
        self._prepare_axonal_collapse()
        self._prepare_dendritic_collapse()

    def _prepare_axonal_collapse(self):
        """Collapse axonal types: IN/PC_DEND_<X>_AX_<Y> → IN/PC_DEND_AX_<Y>"""
        dend_cols = [c for c in self.m_data_df.columns if "_DEND_" in c]
        col_map = {}
        for col in dend_cols:
            match = re.match(r"(IN|PC)_DEND_\d+_AX_(\d+)", col)
            if match:
                cell_type = match.group(1)
                ax_id = match.group(2)
                col_map[col] = f"{cell_type}_DEND_AX_{ax_id}"
        temp_df = self.m_data_df[dend_cols].rename(columns=col_map)
        self.collapsed_axonal_df = temp_df.T.groupby(temp_df.T.index).sum().T

        # Identify excitatory/inhibitory neurons
        self.excitatory_axonal = self.collapsed_axonal_df.columns[
            self.collapsed_axonal_df.columns.str.startswith("PC_DEND_")
        ].tolist()
        self.inhibitory_axonal = self.collapsed_axonal_df.columns[
            self.collapsed_axonal_df.columns.str.startswith("IN_DEND_")
        ].tolist()

    def _prepare_dendritic_collapse(self):
        """Collapse dendritic types: IN/PC_DEND_<X>_AX_<Y> → IN/PC_DEND_<X>"""
        dend_cols = [c for c in self.m_data_df.columns if "_DEND_" in c]
        col_map = {}
        for col in dend_cols:
            match = re.match(r"(IN|PC)_DEND_(\d+)_AX_\d+", col)
            if match:
                cell_type = match.group(1)
                dend_id = match.group(2)
                col_map[col] = f"{cell_type}_DEND_{dend_id}"
        temp_df = self.m_data_df[dend_cols].rename(columns=col_map)
        self.collapsed_dendritic_df = temp_df.T.groupby(temp_df.T.index).sum().T

        # Identify excitatory/inhibitory neurons
        self.excitatory_dendritic = self.collapsed_dendritic_df.columns[
            self.collapsed_dendritic_df.columns.str.startswith("PC_DEND_")
        ].tolist()
        self.inhibitory_dendritic = self.collapsed_dendritic_df.columns[
            self.collapsed_dendritic_df.columns.str.startswith("IN_DEND_")
        ].tolist()

    def plot_layer_distribution(self, df, excitatory, inhibitory, laminar_order, region_name, suffix=""):
        # Only existing layers
        existing_layers = [l for l in laminar_order if l in df.T.columns]
        missing = set(laminar_order) - set(existing_layers)
        if missing:
            print(f"Warning: these layers are missing for {region_name}{suffix}: {missing}")
        df_sel = df.T[existing_layers]

        plt.figure(figsize=(10,6))

        # Excitatory stacked bars
        bottoms = [0]*len(existing_layers)
        for neuron in excitatory:
            if neuron in df_sel.index:
                plt.bar(existing_layers, df_sel.loc[neuron], bottom=bottoms, label=neuron)
                bottoms = [sum(x) for x in zip(bottoms, df_sel.loc[neuron])]

        # Inhibitory lines
        for neuron in inhibitory:
            if neuron in df_sel.index:
                plt.plot(existing_layers, df_sel.loc[neuron], marker='o', linestyle='--', label=neuron)

        plt.xlabel(f"{region_name} Layers")
        plt.ylabel("Neuron density")
        plt.title(f"{region_name} {suffix} collapse: Excitatory (bars) + Inhibitory (lines)")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        filename = os.path.join(self.output_dir, f"{region_name}{suffix}_layer_distribution.png")
        plt.savefig(filename, bbox_inches='tight')
        plt.close()
        print(f"Saved figure: {filename}")

    def generate_all_figures(self):
        regions = {
            "VISp": ["VISp1","VISp23","VISp4","VISp5","VISp6a","VISp6b"],
            "SSs": [col for col in self.m_data_df.T.columns if col.startswith("SSs")],
            "AUDv": [col for col in self.m_data_df.T.columns if col.startswith("AUDv")],
            "ILA": [col for col in self.m_data_df.T.columns if col.startswith("ILA")]
        }

        for region_name, laminar_order in regions.items():
            # First plot axonal collapse
            self.plot_layer_distribution(
                self.collapsed_axonal_df, 
                self.excitatory_axonal, 
                self.inhibitory_axonal, 
                laminar_order, 
                region_name,
                suffix="_axonal"
            )

            # Then plot dendritic collapse
            self.plot_layer_distribution(
                self.collapsed_dendritic_df, 
                self.excitatory_dendritic, 
                self.inhibitory_dendritic, 
                laminar_order, 
                region_name,
                suffix="_dendritic"
            )


# ===========================
# Usage
# ===========================
if __name__ == "__main__":
    analyzer = NeuronLayerAnalyzer(csv_path="data.csv", output_dir="results")
    analyzer.generate_all_figures()
