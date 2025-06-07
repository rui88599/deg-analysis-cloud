
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt

st.set_page_config(page_title="DEG Viewer (Cloud Version)", layout="centered")
st.title("DEG Analysis Viewer (No R required)")

# Upload result CSV (DESeq2 output) and original count matrix
st.subheader("Upload Processed Data")
results_file = st.file_uploader("Upload DESeq2 result CSV (must include log2FoldChange, padj, gene)", type="csv")
counts_file = st.file_uploader("Upload original count matrix (CSV, genes x samples)", type="csv")

padj_cutoff = st.slider("Adjusted p-value cutoff:", 0.01, 1.0, 0.05, 0.01)
logfc_cutoff = st.slider("Log2 fold change cutoff:", 0.0, 5.0, 1.0, 0.1)

if results_file and counts_file:
    result_df = pd.read_csv(results_file)
    counts_df = pd.read_csv(counts_file, index_col=0)

    st.subheader("All DESeq2 Results")
    st.dataframe(result_df.head())

    # Filtering
    filtered_df = result_df[
        (result_df['padj'] < padj_cutoff) &
        (abs(result_df['log2FoldChange']) > logfc_cutoff)
    ].copy()

    st.subheader("Filtered Differentially Expressed Genes")
    st.write(f"Number of genes after filtering: {filtered_df.shape[0]}")
    st.dataframe(filtered_df)

    st.download_button(
        label="Download Filtered Results (CSV)",
        data=filtered_df.to_csv(index=False),
        file_name="filtered_DEGs.csv",
        mime="text/csv"
    )

    # Volcano plot
    result_df["-log10_padj"] = -np.log10(result_df["padj"] + 1e-300)
    result_df["Significance"] = "Not Significant"
    result_df.loc[(result_df["padj"] < padj_cutoff) & (result_df["log2FoldChange"] > logfc_cutoff), "Significance"] = "Upregulated"
    result_df.loc[(result_df["padj"] < padj_cutoff) & (result_df["log2FoldChange"] < -logfc_cutoff), "Significance"] = "Downregulated"

    st.subheader("Volcano Plot")
    fig = px.scatter(
        result_df,
        x="log2FoldChange",
        y="-log10_padj",
        color="Significance",
        hover_name="gene",
        title="Volcano Plot",
        labels={"log2FoldChange": "Log2 Fold Change", "-log10_padj": "-Log10 Adjusted P-value"}
    )
    st.plotly_chart(fig, use_container_width=True)

    # Heatmap
    st.subheader("Heatmap of Top 20 DEGs")
    try:
        top20_genes = filtered_df.sort_values("padj").head(20)["gene"].values
        heatmap_data = counts_df.loc[top20_genes]
        heatmap_data_z = (heatmap_data - heatmap_data.mean(axis=1).values[:, None]) / heatmap_data.std(axis=1).values[:, None]

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.heatmap(heatmap_data_z, cmap="vlag", xticklabels=True, yticklabels=True, ax=ax)
        ax.set_title("Top 20 Differentially Expressed Genes (Z-score normalized)")
        st.pyplot(fig)
    except Exception as e:
        st.warning(f"Heatmap not available: {e}")
else:
    st.info("Please upload both DESeq2 result CSV and count matrix CSV to proceed.")
