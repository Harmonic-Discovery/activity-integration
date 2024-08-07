# Overview

This repo contains code supporting the paper [Leveraging multiple data types for improved compound-kinase bioactivity prediction](https://www.biorxiv.org/content/biorxiv/early/2024/04/16/2024.03.07.583951.full.pdf). 

## Data 

Relevant data is formatted and included in the `data/` directory.

## Single dose -> IC50 prediction with random forests

This repo includes code to perform IC50 prediction from point of concentration data using random forests.  

The script `single_dose_rf.py` performs a simple cross validation search on hyperparameters, and then trains the model on the full training set. 

The script relies on a few standard dependencies, which can be installed via `pip`: `pandas`, `numpy` and `scikit-learn`.

To run the script, use

```bash
python single_dose_rf.py --config single_dose_rf_config.json
```

The config file simply contains the path to the prepped training data and a run name to indentify the model. After training, the final model is saved to a `.pkl` file, which can then be loaded using the `pickle` module.

## Training/prediction with the PWKRR model

Our implementation of pairwise kernel ridge regression is available [here](https://github.com/Harmonic-Discovery/pwkrr); follow the instructions there to install the `pwkrr` package. Model configuration is specified by a `config.json` file (see `pwkrr_config.json` for an example), which allows the user to modulate single vs multi-stage, which split to use, hyperparameters, etc.

The model can be trained using 

```bash
python train_pwkrr.py --config pwkrr_config.json
```

## Software and package versions
- python=3.9
- cython
- pip=22.3.1
- mkl=2021
- blas=1.0=mkl

pip:
- pandas==2.2.2
- rdkit==2023.9.5
- scikit-learn==1.4.2
- scikit-bio==0.6.0
- scipy==1.13.0
- numpy==1.26.4
- joblib==1.4.0


## Reproducing figures

See `figures/` for code/data to reproduce figures appearing in the paper.

## Citation

Pre-print:

```bibtex
@article {Leveraging-multiple-data-types,
	author = {Theisen, Ryan and Wang, Tianduanyi and Ravikumar, Balaguru and Rahman, Rayees and Cicho{\'n}ska, Anna},
	title = {Leveraging multiple data types for improved compound-kinase bioactivity prediction},
	elocation-id = {2024.03.07.583951},
	year = {2024},
	doi = {10.1101/2024.03.07.583951},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/04/16/2024.03.07.583951},
	eprint = {https://www.biorxiv.org/content/early/2024/04/16/2024.03.07.583951.full.pdf},
	journal = {bioRxiv}
}
```
