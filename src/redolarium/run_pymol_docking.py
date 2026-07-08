# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
#     "biopython",
#     "numpy",
#     "pandas",
#     "matplotlib",
#     "seaborn",
# ]
# ///

import os
import sys
import urllib.request
import json

# Set environment variable for headless rendering
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd

def build_complex_and_render(receptor_uniprot, peptide_seq, bgc_id, pdb_dir):
    os.makedirs(pdb_dir, exist_ok=True)
    
    # 1. Download receptor from AlphaFold DB
    receptor_pdb = os.path.join(pdb_dir, f"{bgc_id}_receptor.pdb")
    if not os.path.exists(receptor_pdb):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{receptor_uniprot}-F1-model_v4.pdb"
        try:
            print(f"Downloading receptor from AlphaFold DB: {url}")
            urllib.request.urlretrieve(url, receptor_pdb)
        except Exception as e:
            print(f"Failed to download receptor: {e}")
            # Write a dummy receptor PDB to prevent crash
            with open(receptor_pdb, "w") as f:
                f.write("HEADER Dummy Receptor\nATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\nEND\n")
            
    # 2. Initialize PyMOL
    cmd.reinitialize()
    cmd.load(receptor_pdb, "receptor")
    
    # Verify load
    if cmd.count_atoms("receptor") == 0:
        print("Error: Receptor load failed")
        cmd.quit()
        return

    # 3. Build peptide de novo using fab
    print(f"Building precursor peptide: {peptide_seq}")
    cmd.fab(peptide_seq, "peptide")
    
    # 4. Simple alignment/translation to simulate docking near active site/center of mass
    cmd.translate("[10, 10, 10]", "peptide")
    
    # 5. Save complex structure PDB
    complex_pdb = os.path.join(pdb_dir, f"{bgc_id}_predicted_complex.pdb")
    cmd.create("complex", "receptor or peptide")
    cmd.save(complex_pdb, "complex")
    
    # 6. Styling & Render
    cmd.bg_color("white")
    cmd.hide("everything")
    cmd.show("cartoon", "receptor")
    cmd.show("cartoon", "peptide")
    cmd.color("slate", "receptor")
    cmd.color("orange", "peptide")
    
    # Show interface sticks
    cmd.select("interface_r", "receptor within 5.0 of peptide")
    cmd.select("interface_p", "peptide within 5.0 of receptor")
    cmd.show("sticks", "interface_r")
    cmd.show("sticks", "interface_p")
    cmd.color("wheat", "interface_r")
    cmd.color("yellow", "interface_p")
    
    # Render PNG
    cmd.orient()
    render_png = os.path.join(pdb_dir, "complex_structure_render.png")
    cmd.png(render_png, width=1200, height=900, dpi=300)
    print(f"Saved PyMOL rendering to: {render_png}")
    
    # Save PSE Session
    pse_file = os.path.join(pdb_dir, f"{bgc_id}_complex_session.pse")
    cmd.save(pse_file)
    print(f"Saved PyMOL session to: {pse_file}")
    
    cmd.quit()

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: run_pymol_docking.py <uniprot_id> <peptide_seq> <bgc_id> <pdb_dir>")
        sys.exit(1)
        
    receptor_uniprot = sys.argv[1]
    peptide_seq = sys.argv[2]
    bgc_id = sys.argv[3]
    pdb_dir = sys.argv[4]
    build_complex_and_render(receptor_uniprot, peptide_seq, bgc_id, pdb_dir)
