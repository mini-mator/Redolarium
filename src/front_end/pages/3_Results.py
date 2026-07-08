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
    page_icon="📊",
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
tabular_dir = os.path.join(out_dir, "tabular_data")
bgc_motifs_dir = os.path.join(out_dir, "bgc_motifs")
trees_dir = os.path.join(out_dir, "phylogeny_trees")

num_bgcs = 0
if os.path.exists(bgc_motifs_dir):
    num_bgcs = len([name for name in os.listdir(bgc_motifs_dir) if os.path.isdir(os.path.join(bgc_motifs_dir, name))])

num_trees = 0
if os.path.exists(trees_dir):
    num_trees = len([f for f in os.listdir(trees_dir) if f.endswith('.png') or f.endswith('.svg')])
    
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
st.header("Download Package")

# Zip the directory
if not st.session_state.zip_path:
    with st.spinner("Compressing output package..."):
        zip_base_path = os.path.join(os.path.dirname(out_dir), f"Redolarium_Results_{os.path.basename(out_dir)}")
        zip_file_path = f"{zip_base_path}.zip"
        import zipfile
        try:
            with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                for root, dirs, files in os.walk(out_dir):
                    # Safely skip problematic temp directories
                    if "mmseqs_tmp" in dirs:
                        dirs.remove("mmseqs_tmp")
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
            type="primary",
            use_container_width=True
        )
else:
    st.error("Failed to generate the ZIP archive.")
