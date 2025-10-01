# 🧬 Cancer Subtype Classification from Genomic Mutations

This project develops AI models that classify **cancer subtypes** using **genomic mutation data**.  
Instead of treating mutations as tabular features, we explore **two complementary approaches**:

1. **Tokenized Sequence Approach** — Transform each mutation into a token (e.g. `TP53-MISSENSE-273`)  
   and classify long mutation sequences with **BigBird (Transformer)**.
2. **Tabular Transformer Approach** — Encode mutation features directly as tabular embeddings  
   using a **PyTorch Transformer** and compare with **XGBoost** baselines.

These experiments aim to evaluate how **mutation representation** (sequence vs tabular)  
affects classification performance for complex biological data.

---

## 🚀 Key Experiments
| Notebook | Description |
|----------|--------------|
| `cls_bigbird_protein.ipynb` | Tokenized mutation sequence classification using **BigBird** (long sequence transformer) |
| `cls_transformer_tabular.ipynb` | Tabular feature transformer model + **XGBoost** baseline comparison |
| `preprocessing.ipynb` | Generates token sequences from raw mutation data (`create_token.csv`) |
| `prot_seq_preprocessing.ipynb` | Utilities for applying mutations to UniProt protein sequences |

> **Note**: Original datasets (`train.csv`, `test.csv`) are **not included** in this repository.  
> You can adapt the notebooks to your own mutation data following the same preprocessing pipeline.

---

## 🧰 Tech Stack
- **Python** 3.9+
- **PyTorch**, **Hugging Face Transformers**
- **BigBird**, **nn.Transformer**
- **scikit-learn**, **XGBoost**
- **UniProt REST API** (for optional protein sequence retrieval)

---

## ⚙️ Setup
```bash
# 1. Create a virtual environment
python -m venv venv
source venv/bin/activate

# 2. Install dependencies
pip install -r requirements.txt
# ⚠️ Install PyTorch separately for your system
# See: https://pytorch.org/get-started/locally/
```

---

## 🧪 Running Experiments

### 1️⃣ Tokenized Sequence (BigBird)
```bash
jupyter notebook cls_bigbird_protein.ipynb
```
- Input: `create_token.csv` (generated from preprocessing)
- Model: `BigBirdForSequenceClassification`

### 2️⃣ Tabular Transformer / XGBoost
```bash
jupyter notebook cls_transformer_tabular.ipynb
```
- Input: tabular mutation features
- Models: `nn.Transformer`, `XGBoost`

---

## 🧬 Mutation Preprocessing
- `preprocessing.ipynb` converts mutation annotations into token sequences
- Optional: `prot_seq_preprocessing.ipynb` reconstructs mutated protein sequences using UniProt data

---

## 🧠 Research Goal
This repository investigates whether **tokenizing biological mutations** as textual sequences  
can help **Transformer models** learn more informative representations than traditional tabular features.

---

## 🙏 Acknowledgements
- [Hugging Face Transformers](https://huggingface.co/transformers)
- [PyTorch](https://pytorch.org/)
- [UniProt REST API](https://rest.uniprot.org/)
