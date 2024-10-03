# Improved Protein Localization Predictor

This project uses Python, scikit-learn, and Biopython to predict the subcellular localization of proteins based on their amino acid sequences. The sequences are fetched from the NCBI database using Biopython.

## Features

- Fetches protein sequences from the NCBI database
- Uses character n-grams and Naive Bayes for classification
- Provides accuracy and detailed classification report
- Allows for easy training, evaluation, and prediction of new proteins

## Requirements

- Python 3.6+
- Biopython
- scikit-learn
- numpy

You can install the required packages using pip:

```
pip install biopython scikit-learn numpy
```

## Usage

1. Clone this repository or download the `protein_localization_predictor.py` file.

2. Run the script from the command line, providing your email address as an argument:

   ```
   python protein_localization_predictor.py your.email@example.com
   ```

   Replace `your.email@example.com` with your actual email address. NCBI requires you to provide your email address when using their services.

3. The script will train the model using the provided training data, evaluate it on a test set, and then make predictions for new proteins.

## Customizing the Data

You can easily customize the protein IDs and their localizations used for training and testing. In the `main()` function of the script, modify the following lists:

```python
# Training data
train_protein_ids = ["NP_035357.1", "NP_000207.1", "NP_112445.1", "NP_001320610.1", "NP_001027451.1"]
train_localizations = ["cytoplasm", "nucleus", "mitochondrion", "endoplasmic_reticulum", "plasma_membrane"]

# Test data
test_protein_ids = ["NP_001363856.1", "NP_001365016.1"]
test_localizations = ["golgi_apparatus", "lysosome"]

# New proteins for prediction
new_protein_ids = ["NP_001380.1", "NP_001017963.2"]
```

Replace these lists with your own protein IDs and their corresponding localizations.

## Output

The script will print:
1. The accuracy of the model on the test set
2. A detailed classification report
3. Predictions for the new proteins

## Extending the Project

You can easily extend this project by:
- Adding more protein IDs and localizations to the training and test sets
- Implementing additional machine learning models
- Incorporating more features from the protein sequences

## License

This project is open-source and available under the MIT License.

## Contributing

Contributions, issues, and feature requests are welcome. Feel free to check issues page if you want to contribute.

## Contact

If you have any questions or feedback, please open an issue in this repository.
