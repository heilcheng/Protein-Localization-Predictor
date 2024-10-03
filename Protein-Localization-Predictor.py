import sys
from typing import List, Tuple
from Bio import Entrez, SeqIO
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.metrics import accuracy_score, classification_report
import numpy as np

class ProteinLocalizationPredictor:
    def __init__(self, email: str):
        self.email = email
        Entrez.email = email
        self.vectorizer = CountVectorizer(analyzer='char', ngram_range=(1, 3))
        self.clf = MultinomialNB()

    def fetch_protein_sequences(self, protein_ids: List[str]) -> List[str]:
        sequences = []
        for protein_id in protein_ids:
            try:
                with Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text") as handle:
                    seq_record = SeqIO.read(handle, "genbank")
                    sequences.append(str(seq_record.seq))
            except Exception as e:
                print(f"Error fetching sequence for protein ID {protein_id}: {str(e)}")
        return sequences

    def train(self, protein_ids: List[str], localizations: List[str]) -> None:
        sequences = self.fetch_protein_sequences(protein_ids)
        
        if len(sequences) != len(localizations):
            raise ValueError("Number of sequences and localizations must match")
        
        X = self.vectorizer.fit_transform(sequences)
        self.clf.fit(X, localizations)

    def predict(self, protein_ids: List[str]) -> List[str]:
        sequences = self.fetch_protein_sequences(protein_ids)
        X = self.vectorizer.transform(sequences)
        return self.clf.predict(X)

    def evaluate(self, protein_ids: List[str], localizations: List[str]) -> Tuple[float, str]:
        sequences = self.fetch_protein_sequences(protein_ids)
        X = self.vectorizer.transform(sequences)
        predictions = self.clf.predict(X)
        accuracy = accuracy_score(localizations, predictions)
        report = classification_report(localizations, predictions)
        return accuracy, report

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <your_email@example.com>")
        sys.exit(1)

    email = sys.argv[1]
    predictor = ProteinLocalizationPredictor(email)

    # Training data
    train_protein_ids = ["NP_035357.1", "NP_000207.1", "NP_112445.1", "NP_001320610.1", "NP_001027451.1"]
    train_localizations = ["cytoplasm", "nucleus", "mitochondrion", "endoplasmic_reticulum", "plasma_membrane"]

    # Test data
    test_protein_ids = ["NP_001363856.1", "NP_001365016.1"]
    test_localizations = ["golgi_apparatus", "lysosome"]

    try:
        print("Training the model...")
        predictor.train(train_protein_ids, train_localizations)

        print("\nEvaluating the model...")
        accuracy, report = predictor.evaluate(test_protein_ids, test_localizations)
        print(f"Accuracy: {accuracy:.2f}")
        print("Classification Report:")
        print(report)

        print("\nMaking predictions for new proteins...")
        new_protein_ids = ["NP_001380.1", "NP_001017963.2"]
        predictions = predictor.predict(new_protein_ids)
        for protein_id, prediction in zip(new_protein_ids, predictions):
            print(f"Protein {protein_id}: Predicted localization - {prediction}")

    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
