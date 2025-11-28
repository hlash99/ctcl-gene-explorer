import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import scipy.sparse

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="CTCL Clinical Atlas",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- CUSTOM CSS ---
st.markdown("""
    <style>
    .main {
        background-color: #f8f9fa;
    }
    h1 {
        color: #2c3e50;
        font-family: 'Helvetica Neue', sans-serif;
    }
    .stButton>button {
        width: 100%;
        border-radius: 5px;
        height: 3em;
    }
    .stRadio > div {
        flex-direction: row;
    }
    </style>
    """, unsafe_allow_html=True)

# --- 1. DATA LOADING (CACHED) ---
@st.cache_resource
def load_data():
    """Loads the processed Single-Cell data once."""
    # POINTS TO THE LIGHT FILE NOW
    adata = sc.read_h5ad('ctcl_app_light.h5ad')
    return adata

with st.spinner('Loading the CTCL Atlas...'):
    adata = load_data()

# --- 2. SIDEBAR CONTROLS ---
st.sidebar.title("🧬 CTCL Explorer")
st.sidebar.info("Compare **Tumor** vs. **Eczema** vs. **Normal**.")

st.sidebar.header("Gene Selection")

if 'gene' not in st.session_state:
    st.session_state['gene'] = "TOX"

# Quick Select Buttons
st.sidebar.markdown("### ⚡ Quick Select")
col1, col2 = st.sidebar.columns(2)

if col1.button("TOX (Exhaustion)"):
    st.session_state['gene'] = "TOX"
if col2.button("CCR4 (Target)"):
    st.session_state['gene'] = "CCR4"
if col1.button("CD3E (T-Cell)"):
    st.session_state['gene'] = "CD3E"
if col2.button("MKI67 (Prolif.)"):
    st.session_state['gene'] = "MKI67"

gene_input = st.sidebar.text_input(
    "Or type a gene symbol:", 
    value=st.session_state['gene']
).upper()

st.sidebar.markdown("---")
st.sidebar.caption(
    f"**Dataset Stats:**\n"
    f"• Cells: {adata.n_obs}\n"
    f"• Genes: {adata.n_vars} (Filtered)\n"
    f"• Source: GSE128531"
)

# --- 3. MAIN DASHBOARD ---
st.title(f"Gene Expression: *{gene_input}*")

# FIX 1: Check in the MAIN data, not RAW (since raw is deleted)
if gene_input not in adata.var_names:
    st.error(f"❌ Gene '{gene_input}' not found. Note: This lightweight app only contains top ~5,000 variable genes and clinical markers.")
else:
    # --- ROW 1: THE VISUALIZATION (UMAP) ---
    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("📍 Single-Cell Landscape")
        
        view_option = st.radio(
            "Color Cells By:", 
            ["Gene Expression", "Clinical Condition"], 
            horizontal=True
        )
        
        fig_umap, ax = plt.subplots(figsize=(8, 6))
        
        if view_option == "Gene Expression":
            # FIX 2: Removed use_raw=True because we deleted raw
            sc.pl.umap(
                adata, 
                color=gene_input, 
                ax=ax, 
                show=False, 
                frameon=False, 
                cmap='viridis', 
                size=20, 
                title=f"{gene_input} Expression"
            )
        else:
            sc.pl.umap(
                adata, 
                color='condition', 
                ax=ax, 
                show=False, 
                frameon=False, 
                palette='Set1', 
                size=20,
                title="Clinical Diagnosis"
            )
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

        st.pyplot(fig_umap)

    with col2:
        st.subheader("📊 Clinical Comparison")
        st.markdown("Expression levels by Diagnosis.")
        
        # FIX 3: Handle Sparse Matrix + No Raw
        # Since adata.X is now sparse (Swiss cheese), we must use .toarray() before flattening
        try:
            expression_values = adata[:, gene_input].X.toarray().flatten()
        except:
            # Fallback if it's already dense
            expression_values = adata[:, gene_input].X.flatten()

        df_plot = pd.DataFrame({
            'Expression': expression_values,
            'Condition': adata.obs['condition']
        })
        
        fig_vio, ax = plt.subplots(figsize=(5, 6))
        sns.violinplot(
            data=df_plot, 
            x='Condition', 
            y='Expression', 
            ax=ax, 
            palette='Set2',
            linewidth=1
        )
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("")
        sns.despine()
        st.pyplot(fig_vio)

    # --- ROW 2: ANALYSIS LOGIC ---
    st.markdown("---")
    st.subheader("💡 Analysis Interpretation")
    
    try:
        mean_tumor = df_plot[df_plot['Condition'] == 'CTCL (Tumor)']['Expression'].mean()
        mean_eczema = df_plot[df_plot['Condition'] == 'Eczema (Benign)']['Expression'].mean()
        
        if mean_tumor > mean_eczema:
            diff = mean_tumor - mean_eczema
            if diff > 0.5:
                trend = "SIGNIFICANTLY HIGHER"
                color = "#d9534f" # Red
            else:
                trend = "Slightly Higher"
                color = "#f0ad4e" # Orange
        else:
            trend = "Lower or Same"
            color = "#5cb85c" # Green

        st.markdown(f"""
        <div style="padding: 15px; border-left: 5px solid {color}; background-color: #ffffff;">
            <strong>Automated Insight for {gene_input}:</strong>
            <ul>
                <li>Average expression in <strong>CTCL Tumor</strong>: {mean_tumor:.2f}</li>
                <li>Average expression in <strong>Benign Eczema</strong>: {mean_eczema:.2f}</li>
                <li>Clinical Trend: <span style='color:{color}; font-weight:bold'>{trend}</span> in Tumor vs. Benign mimic.</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        
    except Exception as e:
        st.info("Insufficient data for automated insight.")