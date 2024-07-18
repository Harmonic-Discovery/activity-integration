import numpy as np
import argparse

import json
import os

from sklearn.model_selection import train_test_split, KFold, GridSearchCV
from sklearn.ensemble import RandomForestRegressor

import pandas as pd
import pickle

RANDOM_STATE = 42

def save_pkl(data, filepath):
    """
    Save pkl file

    Arguments:
        data (obj): data to save
        filepath (str): filepath for save location
    """
    with open(filepath, "wb") as fp:
        pickle.dump(data, fp, protocol=pickle.HIGHEST_PROTOCOL)

class SDDataset:
    """
    Utility class that is used for ML models that use single dose observations as features
    """

    def __init__(
        self,
        dataset,
        labels_column="activity_value",
        pivot="concentration",
        limit_obs=-1,
    ):
        """
        Arguments:
            dataset (DataFrame): pandas dataframe obtained from the activity_datasets module
            labels_column (str): if a dataset is passed, the column where the labels are stored (if any)
            pivot (str, default='concentration'): the pivot for the single dose observations
                                                    (either concentration or inhibition)
            limit_obs (int, default=-1): limit to observations that have no more than n_bins + limit obs 
                                            missing concentrations

        Attributes:
            data (array): array containing fingerprints
            labels (array): array containing labels
            targets (array): array of uniprot ids
            smiles (array): list containing smiles strings extracted from the dataset
        """
        assert pivot in [
            "concentration",
            "inhibition",
        ], "pivot needs to be either concentration or inhibition"
        if pivot == "concentration":
            cols = [c for c in dataset.columns if "concentration" in c]
        else:
            cols = [c for c in dataset.columns if "inhibition" in c]
        dataset = dataset[dataset[cols].isnull().sum(axis=1) < (len(cols) + limit_obs)]
        self.data = dataset[cols].to_numpy()
        self.data = np.nan_to_num(self.data, nan=-20)

        self.smiles = dataset["canonical_smiles"].to_numpy()

        if labels_column in dataset.columns:
            self.labels = dataset[labels_column].to_numpy()
        else:
            self.labels = np.array([None] * len(self.smiles))
            print("No labels found in dataset")

        self.targets = dataset["uniprot_id"].to_numpy()
        self.smiles = dataset["canonical_smiles"].to_numpy()

    def __len__(self):
        return self.data.shape[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", nargs="?", default="configs/example.json", type=str)
    args = parser.parse_args()

    os.makedirs("checkpoints", exist_ok=True)
    with open(args.config) as fp:
        config = json.load(fp)

    assert config.get("run_name") is not None

    # load data
    dataset = pd.read_csv(config.get("data_path"))

    dataset = SDDataset(
        dataset, labels_column="activity_value", pivot="concentration", limit_obs=-1
    )
    X, y = dataset.data, dataset.labels

    print(f"{X.shape} total number of train+test samples")
    # train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.15, shuffle=True, random_state=RANDOM_STATE
    )

    # parameter grid for cross validating hyperparameters
    param_grid = {
        "n_estimators": np.arange(50, 400, 50),
        "ccp_alpha": np.exp(np.linspace(-5, 1, 6)),
        "max_depth": np.arange(3, 15, 1),
    }

    print("performing gridsearch")
    model_cv = RandomForestRegressor()
    kf = KFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE).split(X_train)
    gridsearch = GridSearchCV(
        model_cv,
        param_grid=param_grid,
        scoring="neg_root_mean_squared_error",  # "r2",
        cv=kf,
    ).fit(X_train, y_train)

    # fit model with best parameters
    print("fitting best model")
    model = RandomForestRegressor(**gridsearch.best_params_).fit(X_train, y_train)

    model_meta = {
        "model": model,
        "config": config,
    }

    save_pkl(model_meta, f"checkpoints/{config.get('run_name')}.pkl")
