#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import matplotlib.pyplot as plt

class NeuronLayerAnalyzer:
    def __init__(self, csv_path="data.csv", output_dir="results"):
        self.csv_path = csv_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        # Load CSV directly
        self.m_data_df = pd.read_csv(self.csv_path, index_col=0)

        # Identify excitatory/inhibitory neurons (directly from columns)
        self.excitatory = self.m_data_df.columns[self.m_data_df.columns.str.startswith("PC_DEND_")].tolist()
        self.inhibitory = self.m_data_df.columns[self.m_data_df.columns.str.startswith("IN_DEND_")].tolist()

    def plot_region(self, laminar_order, region_name):
        existing_layers = [layer for layer in laminar_order if layer in self.m_data_df.T.columns]
        missing = set(laminar_order) - set(existing_layers)
        if missing:
            print(f"Warning: these layers are missing and will be skipped for {region_name}: {missing}")
    
        df = self.m_data_df.T[existing_layers]
    
        plt.figure(figsize=(10,6))
    
        # Excitatory stacked bars
        bottoms = [0]*len(existing_layers)
        for neuron in self.excitatory:
            if neuron in df.index:
                plt.bar(existing_layers, df.loc[neuron], bottom=bottoms, label=neuron)
                bottoms = [sum(x) for x in zip(bottoms, df.loc[neuron])]
    
        # Inhibitory line plots
        for neuron in self.inhibitory:
            if neuron in df.index:
                plt.plot(existing_layers, df.loc[neuron], marker='o', linestyle='--', label=neuron)
    
        plt.xlabel(f"{region_name} Layers")
        plt.ylabel("Neuron density")
        plt.title(f"Excitatory (bars) + Inhibitory (lines) neuron distribution across {region_name}")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
        # Save figure with full content
        filename = os.path.join(self.output_dir, f"{region_name}.png")
        plt.savefig(filename, bbox_inches='tight')  # <--- key change
        plt.close()
        print(f"Saved figure: {filename}")


    def generate_all_figures(self):
        """Generate plots for all standard regions."""
        regions = {
            "VISp": ["VISp1","VISp23","VISp4","VISp5","VISp6a","VISp6b"],
            "SSs": [col for col in self.m_data_df.T.columns if col.startswith("SSs")],
            "AUDv": [col for col in self.m_data_df.T.columns if col.startswith("AUDv")],
            "ILA": [col for col in self.m_data_df.T.columns if col.startswith("ILA")]
        }

        for region_name, laminar_order in regions.items():
            self.plot_region(laminar_order, region_name)


# ===========================
# Usage
# ===========================
if __name__ == "__main__":
    analyzer = NeuronLayerAnalyzer(csv_path="data.csv", output_dir="results")
    analyzer.generate_all_figures()
