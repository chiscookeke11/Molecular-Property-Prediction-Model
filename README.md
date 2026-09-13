# Molecular Property Prediction Model

A compact, notebook-based baseline for predicting **Tox21 toxicity assay outcomes** from molecular structure. The project turns SMILES strings into a small set of RDKit physicochemical descriptors, then trains a separate class-balanced logistic-regression classifier for each of the 12 Tox21 endpoints.

This repository is intended as an approachable starting point for molecular machine learning: it keeps feature engineering and model evaluation explicit, while accounting for the fact that Tox21 labels are missing independently for each assay.

> **Important:** This is an educational baseline, not a validated toxicity-prediction system. Do not use its output for clinical, regulatory, safety-critical, or chemical-development decisions without appropriate domain validation.

## What is included

| Path | Purpose |
| --- | --- |
| `molecular_property_prediction.ipynb` | End-to-end analysis: data inspection, descriptor calculation, model training, and comparison across endpoints. |
| `sample_data/tox21.csv` | Tox21-style input data containing compound identifiers, SMILES strings, and assay labels. |

## Workflow

The notebook follows this sequence:

1. Load the CSV file with pandas.
2. Inspect dimensions, column types, and endpoint-specific missing labels.
3. Parse every `smiles` value with RDKit. Molecules that RDKit cannot parse are excluded from descriptor-based modeling.
4. Calculate seven molecular descriptors for each valid molecule.
5. For each toxicity endpoint, retain only rows with a known label for that endpoint.
6. Create a stratified 80/20 train/test split (`random_state=24`).
7. Standardize descriptor values with `StandardScaler`, fit on the training split.
8. Fit a `LogisticRegression` model with `class_weight="balanced"`, `max_iter=1000`, and `random_state=2`.
9. Report accuracy, ROC-AUC, a classification report, and a confusion matrix; compare every endpoint with five-fold stratified ROC-AUC cross-validation.

## Dataset format

`sample_data/tox21.csv` has 7,831 molecular records and 14 columns:

- `mol_id`: compound identifier.
- `smiles`: molecular structure encoded as a SMILES string.
- Twelve binary toxicity labels. `1` denotes an active/positive assay result, `0` a negative result, and an empty value means that compound was not tested for that specific endpoint.

The target columns are:

| Family | Endpoint | Description |
| --- | --- | --- |
| Nuclear receptor | `NR-AR` | Androgen receptor |
| Nuclear receptor | `NR-AR-LBD` | Androgen receptor ligand-binding domain |
| Nuclear receptor | `NR-AhR` | Aryl hydrocarbon receptor |
| Nuclear receptor | `NR-Aromatase` | Aromatase |
| Nuclear receptor | `NR-ER` | Estrogen receptor |
| Nuclear receptor | `NR-ER-LBD` | Estrogen receptor ligand-binding domain |
| Nuclear receptor | `NR-PPAR-gamma` | Peroxisome proliferator-activated receptor gamma |
| Stress response | `SR-ARE` | Antioxidant response element |
| Stress response | `SR-ATAD5` | ATPase family AAA domain-containing protein 5 |
| Stress response | `SR-HSE` | Heat-shock response element |
| Stress response | `SR-MMP` | Mitochondrial membrane potential |
| Stress response | `SR-p53` | p53 pathway |

### Missing labels

Missing labels are expected in multi-assay screening data. The notebook **does not drop incomplete rows globally**. Instead, it constructs a different non-null mask for each endpoint, so a compound missing an `NR-AR` result can still contribute to, for example, the `SR-p53` model.

## Molecular features

For a valid RDKit molecule, the notebook calculates the following descriptors:

| Feature | RDKit calculation | Intuition |
| --- | --- | --- |
| `MolWt` | `Descriptors.ExactMolWt` | Exact molecular weight |
| `LogP` | `Descriptors.MolLogP` | Estimated octanol/water partition coefficient |
| `TPSA` | `CalcTPSA` | Topological polar surface area |
| `HBD` | `CalcNumHBD` | Hydrogen-bond donor count |
| `HBA` | `CalcNumHBA` | Hydrogen-bond acceptor count |
| `RotatableBonds` | `CalcNumRotatableBonds` | Molecular flexibility proxy |
| `AromaticRings` | `CalcNumAromaticRings` | Aromatic ring count |

These inexpensive, interpretable features are useful for demonstrating the pipeline but cannot capture all structural information. Fingerprints, graph neural networks, and carefully curated assay-specific features are natural next steps.

## Quick start

### 1. Create an environment

Python 3.10+ is recommended. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate              # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install jupyter pandas numpy scikit-learn rdkit
```

If your platform does not provide an RDKit wheel through pip, use conda instead:

```bash
conda create -n tox21 python=3.11 rdkit pandas numpy scikit-learn jupyter -c conda-forge
conda activate tox21
```

### 2. Start Jupyter

```bash
jupyter notebook molecular_property_prediction.ipynb
```

Run the notebook cells in order. It expects to be launched from the repository root so that `./sample_data/tox21.csv` resolves correctly. The notebook includes a `%pip install rdkit --quiet` cell; installing dependencies in the environment first is still recommended for a repeatable setup.

### 3. Read the results

The `train_and_evaluate(target)` function trains one endpoint model and returns a dictionary containing the fitted model, scaler, sample counts, and evaluation metrics. The final cells iterate over all 12 targets and create a leaderboard sorted by held-out test ROC-AUC.

To train only one endpoint after running the setup cells:

```python
result_nr_ar = train_and_evaluate("NR-AR")
print(result_nr_ar["test_auc"])
```

## Recorded baseline results

The notebook's saved output was produced from the bundled data after eight unparsable SMILES rows were removed. The table below reports the displayed held-out test ROC-AUC and the mean five-fold CV ROC-AUC; it is a baseline snapshot, not a benchmark claim.

| Endpoint | Test ROC-AUC | 5-fold CV ROC-AUC (mean ± SD) |
| --- | ---: | ---: |
| `SR-MMP` | 0.827 | 0.840 ± 0.013 |
| `NR-AhR` | 0.799 | 0.820 ± 0.019 |
| `NR-Aromatase` | 0.786 | 0.798 ± 0.021 |
| `SR-p53` | 0.774 | 0.753 ± 0.021 |
| `NR-ER-LBD` | 0.759 | 0.727 ± 0.030 |
| `NR-AR-LBD` | 0.757 | 0.769 ± 0.036 |
| `NR-AR` | 0.737 | 0.754 ± 0.032 |
| `SR-ARE` | 0.723 | 0.701 ± 0.018 |
| `SR-ATAD5` | 0.707 | 0.699 ± 0.013 |
| `SR-HSE` | 0.687 | 0.700 ± 0.013 |
| `NR-PPAR-gamma` | 0.673 | 0.721 ± 0.030 |
| `NR-ER` | 0.642 | 0.674 ± 0.017 |

Because the endpoints are substantially imbalanced, accuracy alone can be misleading. ROC-AUC and the per-class precision/recall report should be the primary diagnostics in this notebook.

## Reproducibility and limitations

- **Randomness:** the holdout split, classifier, and cross-validation splitter use fixed seeds, so reruns using compatible dependency versions should be broadly reproducible.
- **Invalid structures:** RDKit parsing failures are removed before modeling. Inspect and curate these records rather than silently discarding them in a production workflow.
- **Class imbalance:** `class_weight="balanced"` reduces majority-class bias, but it does not replace threshold tuning, calibrated probabilities, or appropriate precision-recall analysis.
- **Cross-validation preprocessing:** the current cross-validation call scales the full endpoint feature matrix before splitting folds. That lets validation-fold distribution information influence scaling. For a rigorously unbiased CV estimate, replace it with an sklearn `Pipeline` that contains `StandardScaler` and `LogisticRegression`, then pass that pipeline to `cross_val_score`.
- **Chemical generalization:** the random split may place close structural analogues in both train and test sets. Use scaffold-aware splits and an external test set to assess real-world generalization.
- **Scope:** seven handcrafted descriptors and a linear classifier are deliberately simple. Performance may improve with molecular fingerprints, hyperparameter tuning, feature selection, or graph-based models, but those changes require careful validation.

## Suggested improvements

1. Use an sklearn `Pipeline` for preprocessing and model fitting during cross-validation.
2. Add Morgan/ECFP fingerprints and compare them with the descriptor baseline.
3. Evaluate precision-recall AUC, balanced accuracy, calibration, and endpoint-specific decision thresholds.
4. Use Bemis-Murcko scaffold splits to reduce structural leakage.
5. Save trained pipelines and descriptor metadata for inference on new SMILES.
6. Add automated tests for SMILES parsing, descriptor column order, and per-endpoint missing-label handling.

## License and data provenance

No license file or standalone dataset provenance notice is currently included in this repository. Before redistributing the data or using the project beyond experimentation, verify the applicable source-data terms and add an explicit project license.
