import streamlit as st
import sys
import os
import shutil
from pathlib import Path

_src_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _src_dir not in sys.path:
    sys.path.insert(0, _src_dir)

from front_end.utils.state_manager import init_session_state

st.set_page_config(
    page_title="Redolarium | Results",
    page_icon="",
    layout="wide"
)

init_session_state()

st.markdown("""
<style>
    .reportview-container {
        background: #f8f9fa;
    }
    h1, h2, h3 {
        color: #0b3d91;
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    }
    .disclaimer {
        font-size: 0.9em;
        color: #a00;
        border-left: 4px solid #a00;
        padding-left: 10px;
        margin-bottom: 20px;
        background-color: #fff0f0;
        padding: 10px;
    }
    .metric-card {
        background: white;
        padding: 15px;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
        border-top: 4px solid #0b3d91;
    }
    .metric-value {
        font-size: 2em;
        font-weight: bold;
        color: #0b3d91;
    }
    .metric-label {
        font-size: 0.9em;
        color: #555;
        text-transform: uppercase;
    }
</style>
""", unsafe_allow_html=True)

st.title("Analysis Results")

st.markdown("""
<div class='disclaimer'>
    <strong>CRITICAL NOTICE - DATA RETENTION POLICY:</strong> All results and submitted sequences are 
    permanently deleted from our servers 7 days after execution. Please download your results immediately.
</div>
""", unsafe_allow_html=True)

if not st.session_state.tmp_out_dir or not os.path.exists(st.session_state.tmp_out_dir):
    st.error("No result directory found. Did the pipeline run successfully?")
    st.stop()

out_dir = st.session_state.tmp_out_dir

# Calculate basic metrics based on output folder
bgc_base_dir = os.path.join(out_dir, "bgc")
bgc_list = []
num_bgcs = 0
num_trees = 0

if os.path.exists(bgc_base_dir):
    bgc_list = sorted([d for d in os.listdir(bgc_base_dir) if os.path.isdir(os.path.join(bgc_base_dir, d))])
    num_bgcs = len(bgc_list)
    for bgc_name in bgc_list:
        trees_dir = os.path.join(bgc_base_dir, bgc_name, "phylogeny_trees")
        if os.path.exists(trees_dir):
            num_trees += len([f for f in os.listdir(trees_dir) if f.endswith('.png') or f.endswith('.svg')])
            
total_files = sum([len(files) for r, d, files in os.walk(out_dir)])
total_size = sum(os.path.getsize(os.path.join(dirpath, filename)) for dirpath, _, filenames in os.walk(out_dir) for filename in filenames)
size_mb = total_size / (1024 * 1024)

st.header("Execution Summary")
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.markdown(f"<div class='metric-card'><div class='metric-value'>{num_bgcs}</div><div class='metric-label'>BGCs Analyzed</div></div>", unsafe_allow_html=True)
with col2:
    st.markdown(f"<div class='metric-card'><div class='metric-value'>{num_trees}</div><div class='metric-label'>Phylogenetic Trees</div></div>", unsafe_allow_html=True)
with col3:
    st.markdown(f"<div class='metric-card'><div class='metric-value'>{total_files}</div><div class='metric-label'>Files Generated</div></div>", unsafe_allow_html=True)
with col4:
    st.markdown(f"<div class='metric-card'><div class='metric-value'>{size_mb:.2f} MB</div><div class='metric-label'>Total Data Volume</div></div>", unsafe_allow_html=True)

st.markdown("---")

if num_bgcs > 0:
    st.header("BGC Interactive Gallery")
    selected_bgc = st.selectbox("Select a BGC to view generated insights:", bgc_list)
    
    if selected_bgc:
        bgc_target_dir = os.path.join(bgc_base_dir, selected_bgc)
        
        tab1, tab2, tab3, tab4 = st.tabs(["Synteny Map", "Metabolic Linkage", "Environmental Cues", "Phylogeny"])
        
        with tab1:
            st.subheader(f"{selected_bgc} - Dual-Scale Synteny Map")
            synteny_img = os.path.join(bgc_target_dir, "hgt_evolution", f"{selected_bgc}_hgt_synteny.png")
            if os.path.exists(synteny_img):
                st.image(synteny_img, use_column_width=True)
            else:
                st.info("Synteny map not found or not generated for this cluster.")
                
        with tab2:
            st.subheader(f"{selected_bgc} - Pathway Linkage Network")
            linkage_img = os.path.join(bgc_target_dir, "metabolic_pathways", f"{selected_bgc}_metabolic_linkage.png")
            if os.path.exists(linkage_img):
                st.image(linkage_img, use_column_width=True)
            else:
                st.info("Metabolic linkage map not found.")
                
        with tab3:
            st.subheader(f"{selected_bgc} - Genome-Aware Cues")
            cues_csv = os.path.join(bgc_target_dir, "regulatory_networks", f"{selected_bgc}_environmental_cues.csv")
            if os.path.exists(cues_csv):
                import pandas as pd
                st.dataframe(pd.read_csv(cues_csv))
            else:
                st.info("No strong environmental cues or promoters mapped for this cluster.")
                
        with tab4:
            st.subheader(f"{selected_bgc} - Clocking & Trees")
            trees_dir = os.path.join(bgc_target_dir, "phylogeny_trees")
            if os.path.exists(trees_dir):
                trees = [f for f in os.listdir(trees_dir) if f.endswith('.png')]
                if trees:
                    for t in trees:
                        st.image(os.path.join(trees_dir, t), caption=t, use_column_width=True)
                else:
                    st.info("No phylogenetic trees available.")
            else:
                st.info("No phylogenetic trees available.")

st.markdown("---")
st.markdown("---")
st.header("Download Package")

# Zip the directory
zip_base_path = os.path.join(os.path.dirname(out_dir), f"Redolarium_Results_{os.path.basename(out_dir)}")
zip_file_path = f"{zip_base_path}.zip"

if not st.session_state.zip_path or not os.path.exists(zip_file_path):
    with st.spinner("Compressing output package..."):
        import zipfile
        try:
            with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                for root, dirs, files in os.walk(out_dir):
                    # Safely skip problematic temp directories
                    if "mmseqs_tmp" in dirs:
                        dirs.remove("mmseqs_tmp")
                    if "antismash_results" in dirs:
                        dirs.remove("antismash_results")
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.relpath(file_path, out_dir)
                        try:
                            zf.write(file_path, arcname)
                        except OSError:
                            pass # Skip files that are locked by the system
            st.session_state.zip_path = zip_file_path
        except Exception as e:
            st.error(f"Failed to zip results: {e}")
            st.session_state.zip_path = None

if st.session_state.zip_path and os.path.exists(st.session_state.zip_path):
    with open(st.session_state.zip_path, "rb") as fp:
        btn = st.download_button(
            label="Download Full Results (.zip)",
            data=fp,
            file_name=os.path.basename(st.session_state.zip_path),
            mime="application/zip",
            type="primary"

        )
else:
    st.error("Failed to generate the ZIP archive.")
