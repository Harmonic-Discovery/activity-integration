from pwkrr import PairwiseKernelRidge, PairwiseKernelRidgeDataset
from scipy.stats import pearsonr, spearmanr, kendalltau

import os
import time
import pandas as pd
import numpy as np
import json
import argparse

RANDOM_STATE = 42

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--config", nargs="?", default="example_config.json", type=str)
    parser.add_argument("--test", dest="test", action="store_true")
    parser.set_defaults(test=False)
    args = parser.parse_args()

    os.makedirs("checkpoints", exist_ok=True)
    with open(args.config) as fp:
        params = json.load(fp)

    assert "run_name" in params.keys()
    os.makedirs(os.path.join("checkpoints", params.get("run_name")), exist_ok=True)

    with open(
        os.path.join("checkpoints", params.get("run_name"), "params.json"), "w"
    ) as fp:
        json.dump(params, fp)

    split = params.get("split", "ck")
    stage = params.get("stage", "single")
    dataset = pd.read_csv(params.get("data_path"))
    dataset = dataset.dropna()

    if split == "ck":
        val_dataset = dataset[dataset["ck_split_set"] == "test"]
        train_dataset = dataset[dataset["ck_split_set"] == "train"]
        if stage == "single":
            train_dataset = train_dataset[train_dataset["activity_type"] == "measured"]
            runid = "pwkrr_ck_single"
        else:
            runid = "pwkrr_ck_multi"
    elif split == "ligand":
        val_dataset = dataset[dataset["ligand_split_set"] == "test"]
        train_dataset = dataset[dataset["ligand_split_set"] == "train"]
        if stage == "single":
            train_dataset = train_dataset[train_dataset["activity_type"] == "measured"]
            runid = "pwkrr_ligand_single"
        else:
            runid = "pwkrr_ligand_multi"
    elif split == "scaffold":
        val_dataset = dataset[dataset["scaffold_split_set"] == "test"]
        train_dataset = dataset[dataset["scaffold_split_set"] == "train"]
        if stage == "single":
            train_dataset = train_dataset[train_dataset["activity_type"] == "measured"]
            runid = "pwkrr_scaffold_single"
        else:
            runid = "pwkrr_scaffold_multi"
    elif split == "cluster":
        val_dataset = dataset[dataset["cluster_split_set"] == "test"]
        train_dataset = dataset[dataset["cluster_split_set"] == "train"]
        if stage == "single":
            train_dataset = train_dataset[train_dataset["activity_type"] == "measured"]
            runid = "pwkrr_cluster_single"
        else:
            runid = "pwkrr_cluster_multi"
    elif split == "prod":
        val_dataset = None
        train_dataset = dataset
        if stage == "single":
            train_dataset = train_dataset[train_dataset["activity_type"] == "measured"]
            runid = "pwkrr_prod_single"
        else:
            runid = "pwkrr_prod_multi"

    if args.test:
        # for the sake of testing, train a model on a tiny subset
        train_dataset = train_dataset.sample(n=100)
        val_dataset = val_dataset.sample(n=100)

    train_dataset = PairwiseKernelRidgeDataset(
        train_dataset, labels_column="activity_value", n_workers=4, params=params
    )

    if val_dataset is not None:
        val_dataset = PairwiseKernelRidgeDataset(
            val_dataset, labels_column="activity_value", n_workers=4, params=params
        )

    model = PairwiseKernelRidge(
        alpha=params.get("alpha", 0.5),
        protein_kernel=params.get("protein_kernel", "normalized_ssw"),
        ligand_kernel=params.get("ligand_kernel", "tanimoto"),
        protein_kernel_power=params.get("protein_kernel_power", 1.0),
        ligand_kernel_power=params.get("ligand_kernel_power", 1.0),
    )
    tic = time.time()
    print("Training...")
    model = model.fit(
        proteins=train_dataset.proteins,
        ligands=train_dataset.ligands,
        y=train_dataset.labels,
        maxiter=params.get("max_iter", 1000),
        store_eigs=params.get("store_eigs", False),
    )
    print("Training completed")
    train_time = time.time() - tic

    # EVALUATION
    if val_dataset is not None:
        tic = time.time()
        y_hat = model.predict(val_dataset.proteins, val_dataset.ligands)
        test_time = time.time() - tic

        y_val = val_dataset.labels

        val_rmse = np.sqrt(np.mean((y_hat - y_val) ** 2))
        val_pearson, _ = pearsonr(y_hat, y_val)
        val_spearman, _ = spearmanr(y_hat, y_val)
        val_kendall, _ = kendalltau(y_hat, y_val)

        print("Validation completed")
        print(f"RMSE: {round(val_rmse, 3)}")
        print(f"Pearson: {round(val_pearson, 3)}")
        print(f"Spearman: {round(val_spearman, 3)}")

        meta = {
            "val_rmse": val_rmse,
            "val_pearson": val_pearson,
            "val_spearman": val_spearman,
            "val_kendall": val_kendall,
            "params": params,
        }

    else:
        meta = {"params": params}

    model_meta = model.save(
        os.path.join("checkpoints", params.get("run_name"), "model.pkl"), meta=meta
    )
