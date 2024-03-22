# Protein Localization Predictor

This project uses Python and scikit-learn to predict the subcellular localization of proteins based on their amino acid sequences. The sequences are fetched from the NCBI database using Biopython.

## Requirements

- Python 3
- Biopython
- scikit-learn

## Usage

1. Replace `"your.email@example.com"` with your actual email address in `ProteinLocalizationPredictor.py`. NCBI requires you to provide your email address when using their services.

2. Replace the protein IDs and their localizations in `ProteinLocalizationPredictor.py` with your data. The current protein IDs in the script are `"NP_035357.1"`, `"NP_000207.1"`, and `"NP_112445.1"`, and their corresponding localizations are `"cytoplasm"`, `"nucleus"`, and `"mitochondrion"`.

3. Run `ProteinLocalizationPredictor.py`:

```bash
python3 ProteinLocalizationPredictor.py
