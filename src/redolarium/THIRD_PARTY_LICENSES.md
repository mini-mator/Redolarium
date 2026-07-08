# Third-Party Licences — Redolarium v3.0.0

Redolarium is MIT-licensed software that relies on the following third-party tools,
libraries, and databases. Each retains its own licence terms as described below.

---

## Python Libraries

| Package | Version | Licence | Link |
|---------|---------|---------|------|
| Biopython | ≥1.80 | Biopython Licence (BSD-like) | https://biopython.org/wiki/License |
| pyhmmer | ≥0.10.2 | MIT | https://github.com/althonos/pyhmmer |
| pyrodigal | ≥3.0.0 | GPL-3.0 | https://github.com/althonos/pyrodigal |
| numpy | any | BSD-3-Clause | https://numpy.org |
| pandas | any | BSD-3-Clause | https://pandas.pydata.org |
| matplotlib | any | PSF (BSD-like) | https://matplotlib.org |
| seaborn | any | BSD-3-Clause | https://seaborn.pydata.org |
| openpyxl | any | MIT | https://openpyxl.readthedocs.io |
| networkx | any | BSD-3-Clause | https://networkx.org |
| scipy | any | BSD-3-Clause | https://scipy.org |
| streamlit | any | Apache-2.0 | https://streamlit.io |
| python-docx | any | MIT | https://python-docx.readthedocs.io |
| pydantic | any | MIT | https://docs.pydantic.dev |

> **pyrodigal note:** pyrodigal wraps Prodigal (GPL-3.0). This means that
> redistribution of a compiled binary that statically links pyrodigal must comply
> with GPL-3.0. Redolarium itself remains MIT; pyrodigal is a runtime dependency.

---

## External Tools (Optional / Invoked via Subprocess)

### antiSMASH
- **Licence:** GNU Affero General Public License v3.0 (AGPL-3.0)
- **Source:** https://github.com/antismash/antismash
- **Usage in Redolarium:** Optional; invoked via Docker container as primary BGC
  detection engine if Docker/WSL is present. Falls back to internal HMMER-based
  scanner if antiSMASH is unavailable.
- **Note:** antiSMASH is NOT bundled with Redolarium. Users must install it
  separately via the official Docker image (`antismash/standalone`).
- **Citation:** Blin et al. 2023, Nucleic Acids Research (doi:10.1093/nar/gkad344)

### Smina / AutoDock Vina
- **Licence:** Apache License 2.0
- **Source:** https://github.com/mwojcikowski/smina
- **Usage in Redolarium:** Optional; invoked for molecular docking if installed.
- **Citation:** Koes et al. 2013, Journal of Chemical Information and Modeling
  (doi:10.1021/ci300604z)

### PyMOL
- **Licence:** PyMOL Open Source (BSD-like) or Schrödinger commercial
- **Source:** https://github.com/schrodinger/pymol-open-source
- **Usage in Redolarium:** Optional; invoked headlessly for structure rendering.
- **Note:** Only the open-source PyMOL build is expected. Users must install separately.

### Bakta / Prokka (optional annotation engines)
- **Bakta Licence:** GPL-3.0 — https://github.com/oschwengers/bakta
- **Prokka Licence:** GPL-3.0 — https://github.com/tseemann/prokka
- **Usage:** Optional; invoked via Snakemake workflow rules if installed in PATH.

---

## Databases

### MIBiG — Minimum Information about a Biosynthetic Gene cluster
- **Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Version used:** MIBiG 3.1
- **Source:** https://mibig.secondarymetabolites.org/
- **Attribution (required by CC BY 4.0):**
  Terlouw et al. 2023, Nucleic Acids Research (doi:10.1093/nar/gkac1049)
- **Usage:** BGC profile gene lists and compound annotations are embedded in
  `bgc_analysis.py` and referenced during BGC similarity scoring.

### Pfam — Protein Families Database
- **Licence:** Creative Commons Zero (CC0 / Public Domain)
- **Version used:** Pfam-A 36.0 (curated essential subset)
- **Source:** https://www.ebi.ac.uk/interpro/download/pfam/
- **Attribution:**
  Mistry et al. 2021, Nucleic Acids Research (doi:10.1093/nar/gkaa913)
- **Usage:** HMM profiles bundled in `resources/essential_bgc.hmm` and
  `resources/curated_screening.hmm` for BGC and cargo detection.

### NCBI / GenBank / BLAST
- **Licence:** Public domain (US Government work); ToS apply to API usage
- **API ToS:** https://www.ncbi.nlm.nih.gov/home/about/policies/
- **Usage:** Genome download via Entrez eFetch; optional remote BLASTp.
- **Requirement:** A valid user email address is required by NCBI ToS.

### ChEMBL
- **Licence:** Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)
- **Source:** https://www.ebi.ac.uk/chembl/
- **Attribution (required by CC BY-SA 4.0):**
  Mendez et al. 2019, Nucleic Acids Research (doi:10.1093/nar/gky1075)
- **Usage:** Bioactivity (Ki) data fetched at runtime via ChEMBL REST API.
  Results are cached locally and must NOT be redistributed in bundled datasets
  without the CC BY-SA 4.0 licence notice.

### AlphaFold Database
- **Licence:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Source:** https://alphafold.ebi.ac.uk/
- **Attribution (required by CC BY 4.0):**
  Varadi et al. 2022, Nucleic Acids Research (doi:10.1093/nar/gkab1061)
- **Usage:** Receptor structure coordinates fetched at runtime for docking.

### KEGG — Kyoto Encyclopedia of Genes and Genomes
- **Licence:** KEGG licence (free for academic/non-commercial use via REST API;
  commercial use requires subscription)
- **Source:** https://www.genome.jp/kegg/
- **Attribution:** Kanehisa & Goto 2000, Nucleic Acids Research (doi:10.1093/nar/28.1.27)
- **Usage:** Metabolic pathway and EC number lookups via KEGG REST API.
- **IMPORTANT:** Use of the KEGG API for commercial purposes requires a paid
  KEGG FTP/API subscription. Academic researchers may use the free REST API.

---

## Algorithms / Methods

| Method | Reference | Used in |
|--------|-----------|---------|
| Bacterial molecular clock | Ochman et al. 2000 (doi:10.1016/S0092-8674(00)80405-8) | `evolution.py`, `phylogeny.py` |
| antiSMASH BGC clustering gap (10 kb) | Blin et al. 2023 (doi:10.1093/nar/gkad344) | `bgc_analysis.py` |
| Promoter sigma box consensus | Hawley & McClure 1983; Helmann 1995 | `bgc_analysis.py` |
| Shine-Dalgarno detection | Shine & Dalgarno 1974; Steitz & Jakes 1975 | `bgc_analysis.py` |
| Radical SAM sactipeptide (PF04738) | Grell et al. 2018 (doi:10.1038/s41589-018-0122-9) | `bgc_analysis.py` |
| Difficidin EF-Tu mechanism | Broderick et al. 2021 (doi:10.1073/pnas.2019378118) | `bgc_analysis.py` |
| 7-gene MLSA scheme | Goris et al. 2007 (doi:10.1099/ijs.0.65060-0) | `annotation.py` |
| AlphaFold-Multimer | Evans et al. 2021 (doi:10.1101/2021.10.04.463034) | `docking.py` |
