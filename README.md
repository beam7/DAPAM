# DAPAM Database of Antibacterial Peptides and their Mechanisms
[![License: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/Code-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Antibacterial peptides offer promising alternatives to conventional antibiotics, yet no comprehensive dataset previously linked ABP sequences to their mechanisms of action. DAPAM fills this gap by pairing sequences with mechanism annotations including:

- **Pore formation** (toroidal pore, barrel-stave)
- **Carpet model**
- **Membrane disruption**
- **Immunomodulatory effects**
- **Intracellular targeting**


## Pipeline

![DAPAM Pipeline](https://github.com/beam7/DAPAM/blob/main/DAPAMpipeline.png)

### Key Statistics

- **1,716** articles retrieved from PubMed (700 full-text)
- **37.7%** of articles matched the sequence regex
- **324** amino acid sequences extracted
- **283** confirmed antibacterial peptides
- **208** unique ABPs with mechanism annotations
- **70** sequences (34%) unique to DAPAM

## Dataset Files

| File | Description |
|------|-------------|
| `DAPAM.json` | Full dataset with all annotations |
| `DAPAM.fasta` | Sequences in FASTA format |
| `DAPAM.csv` | Tabular format for analysis |

### Data Fields

Each entry contains:

- `sequence`: Amino acid sequence (10-60 residues)
- `mechanism`: Primary mechanism of action
- `mechanism_keywords`: Keywords from source literature
- `bacterial_target`: Target bacterial type (Gram+/Gram-/both)
- `source_sentences`: Evidence sentences from literature
- `source_pmid`: PubMed ID of source article

## Installation

```bash
git clone https://github.com/beam7/DAPAM.git
cd DAPAM
pip install -r requirements.txt
```

### Dependencies

- Python 3.8+
- pandas
- numpy
- scipy
- biopython
- matplotlib
- seaborn
- peptides

## Mechanism Distribution

| Mechanism | Count |
|-----------|-------|
| Pore | 159 |
| Carpet | 40 |
| Membrane disruption | 5 |
| Immunomodulatory | 3 |
| Intracellular targeting | 1 |

## Statistical Validation

DAPAM was compared against **16,439 ABPs from Peptipedia** to confirm representativeness. Results show:

- Most physicochemical features have similar distributions
- Only percent cysteine exceeded the small effect size threshold (Cohen's d > 0.2)
- Literature mining did not introduce significant bias

## Citation

If you use DAPAM in your research, please cite this repo. Currently working to publish.

## License

- **Dataset**: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- **Code**: [MIT License](LICENSE)

## Acknowledgments

B.M. gratefully acknowledges the Research Advancements in Diversity Scholars (RADS) program at Cedars-Sinai Health Language Processing Center.

## Contact

- **Beatrice Mihalache** — [beatrice.mihalache@ucsf.edu](mailto:beatrice.mihalache@ucsf.edu)
- Department of Biochemistry and Biophysics, University of California San Francisco

## Related Resources

- [APD3: Antimicrobial Peptide Database](https://aps.unmc.edu/)
- [DBAASP](https://dbaasp.org/)
- [Peptipedia](https://peptipedia.cl/)
