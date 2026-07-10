# Redolarium
**A Universal Genomic, Metabolic, and Biosynthetic Gene Cluster (BGC) Analysis Pipeline**

Redolarium is an automated bioinformatics pipeline that evaluates sequence similarity, stoichiometric precursors, genome-aware promoter prediction, and horizontal gene transfer (HGT) signatures. It features a robust Python backend with a streamlined Streamlit graphical interface.

## System Requirements
- **OS:** Windows, macOS, or Linux
- **Python:** 3.10+
- **Dependencies:** Listed in `requirements.txt` or `environment.yml`

## Installation
Clone the repository and install the dependencies using pip:
```bash
git clone https://github.com/mini-mator/Redolarium.git
cd Redolarium
pip install -r requirements.txt
```

## Quick Start
To launch the interactive graphical user interface:
```bash
# On Windows, simply double click run.bat, or:
streamlit run src/front_end/app.py
```
Upload your `.gb` sequence (like the one provided in `example_data/`) or provide an NCBI Accession ID to begin the analysis.
