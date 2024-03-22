from Bio import Entrez, SeqIO
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.metrics import accuracy_score

Entrez.email = "your.email@example.com"  #Replace with your email

def fetch_protein_sequences(protein_ids):
    sequences = []
    for protein_id in protein_ids:
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()
        sequences.append(str(seq_record.seq))
    return sequences

protein_ids = ["NP_035357.1", "NP_000207.1", "NP_112445.1"]
localizations = ["cytoplasm", "nucleus", "mitochondrion"] #Replace with protein IDs and localizations

sequences = fetch_protein_sequences(protein_ids)

sequences_train, sequences_test, localizations_train, localizations_test = train_test_split(sequences, localizations, test_size=0.2, random_state=42)

vectorizer = CountVectorizer(analyzer='char', ngram_range=(1, 1))
X_train = vectorizer.fit_transform(sequences_train)
X_test = vectorizer.transform(sequences_test)
clf = MultinomialNB()
clf.fit(X_train, localizations_train)

localizations_pred = clf.predict(X_test)

print("Accuracy:", accuracy_score(localizations_test, localizations_pred))