#!/usr/bin/env python3

import argparse
import gc
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
)


BASE_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

NUMERIC_COLS = [
    "DP",
    "ALT_COUNT",
    "REF_COUNT",
    "VAF",
    "ALT_FWD",
    "ALT_REV",
    "SB",
    "ALT_BQ_MEAN",
    "HOMOPOLY",
    "GC11",
]

FEATURE_NAMES = [
    "DP",
    "ALT_COUNT",
    "REF_COUNT",
    "VAF",
    "ALT_FWD",
    "ALT_REV",
    "SB",
    "ALT_BQ_MEAN",
    "HOMOPOLY",
    "GC11",
    "CTX5_0",
    "CTX5_1",
    "CTX5_2",
    "CTX5_3",
    "CTX5_4",
]

REQUIRED_BASE_COLS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "LABEL",
    "DP",
    "ALT_COUNT",
    "REF_COUNT",
    "VAF",
    "ALT_FWD",
    "ALT_REV",
    "SB",
    "ALT_BQ_MEAN",
    "CTX5",
    "HOMOPOLY",
    "GC11",
]


def parse_sample_file(value):
    """
    Parse SAMPLE=/path/to/features.tsv.gz.
    """
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"Expected SAMPLE=PATH format, got: {value}"
        )

    sample, path = value.split("=", 1)

    sample = sample.strip()
    path = path.strip()

    if not sample:
        raise argparse.ArgumentTypeError(f"Missing sample name in: {value}")

    if not path:
        raise argparse.ArgumentTypeError(f"Missing path in: {value}")

    return sample, path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train SNP and INDEL XGBoost rescue classifiers from feature TSV files."
    )

    parser.add_argument(
        "--train",
        action="append",
        type=parse_sample_file,
        required=True,
        help="Training feature file in SAMPLE=PATH format. Can be used multiple times.",
    )
    parser.add_argument(
        "--valid",
        action="append",
        type=parse_sample_file,
        required=True,
        help="Validation feature file in SAMPLE=PATH format. Can be used multiple times.",
    )
    parser.add_argument(
        "--test",
        action="append",
        type=parse_sample_file,
        required=True,
        help="Test feature file in SAMPLE=PATH format. Can be used multiple times.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where models, metrics, predictions, and threshold tables will be saved.",
    )

    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--nthread", type=int, default=None)

    parser.add_argument("--num-boost-round", type=int, default=2000)
    parser.add_argument("--early-stopping-rounds", type=int, default=75)
    parser.add_argument("--verbose-eval", type=int, default=50)

    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--eta", type=float, default=0.03)
    parser.add_argument("--subsample", type=float, default=0.85)
    parser.add_argument("--colsample-bytree", type=float, default=0.85)
    parser.add_argument("--min-child-weight", type=float, default=5.0)
    parser.add_argument("--reg-lambda", type=float, default=2.0)
    parser.add_argument("--reg-alpha", type=float, default=0.1)
    parser.add_argument("--tree-method", default="hist")

    parser.add_argument(
        "--models",
        choices=["both", "snp", "indel"],
        default="both",
        help="Which model(s) to train [default: both]",
    )

    return parser.parse_args()


def sample_file_list_to_dict(items):
    out = {}

    for sample, path in items:
        if sample in out:
            raise ValueError(f"Duplicate sample name: {sample}")
        out[sample] = path

    return out


def check_input_files(file_map, group_name):
    for sample, path in file_map.items():
        if not Path(path).is_file():
            raise FileNotFoundError(
                f"Missing {group_name} feature file for {sample}: {path}"
            )


def encode_ctx5_array(series):
    arr = []

    for value in series.fillna("NNNNN").astype(str):
        ctx = value.strip().upper()

        if len(ctx) != 5:
            ctx = "NNNNN"

        arr.append([BASE_MAP.get(base, 4) for base in ctx])

    return np.asarray(arr, dtype=np.float32)


def normalize_columns(df, sample, path):
    """
    Supports both:
      - feature files with TYPE column
      - feature files with VAR_TYPE column
      - feature files with or without SAMPLE column
    """

    if "VAR_TYPE" not in df.columns and "TYPE" in df.columns:
        df = df.rename(columns={"TYPE": "VAR_TYPE"})

    if "SAMPLE" not in df.columns:
        df["SAMPLE"] = sample
    else:
        df["SAMPLE"] = sample

    missing = [col for col in REQUIRED_BASE_COLS if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing expected columns: {missing}")

    if "VAR_TYPE" not in df.columns:
        raise ValueError(
            f"{path} is missing variant type column. Expected TYPE or VAR_TYPE."
        )

    keep_cols = [
        "SAMPLE",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "VAR_TYPE",
        "LABEL",
        *NUMERIC_COLS,
        "CTX5",
    ]

    return df[keep_cols].copy()


def load_one(path, sample, want_types):
    print(f"Loading {sample}: {path}")

    df = pd.read_csv(path, sep="\t", compression="gzip")
    df = normalize_columns(df, sample, path)

    df["VAR_TYPE"] = df["VAR_TYPE"].astype(str).str.upper()
    df = df[df["VAR_TYPE"].isin(want_types)].copy()

    df["LABEL"] = pd.to_numeric(df["LABEL"], errors="coerce").fillna(0).astype(np.int8)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").fillna(0).astype(np.int64)

    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def load_group(file_map, want_types):
    frames = []

    for sample, path in file_map.items():
        df = load_one(path, sample, want_types)
        frames.append(df)

    if not frames:
        raise ValueError("No feature files were provided.")

    return pd.concat(frames, ignore_index=True)


def make_matrix(df):
    df = df.copy()

    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(np.float32)

    x = pd.DataFrame(
        {
            "DP": df["DP"].astype(np.float32),
            "ALT_COUNT": df["ALT_COUNT"].astype(np.float32),
            "REF_COUNT": df["REF_COUNT"].astype(np.float32),
            "VAF": df["VAF"].astype(np.float32),
            "ALT_FWD": df["ALT_FWD"].astype(np.float32),
            "ALT_REV": df["ALT_REV"].astype(np.float32),
            "SB": df["SB"].astype(np.float32),
            "ALT_BQ_MEAN": df["ALT_BQ_MEAN"].astype(np.float32),
            "HOMOPOLY": df["HOMOPOLY"].astype(np.float32),
            "GC11": df["GC11"].astype(np.float32),
        }
    )

    ctx = encode_ctx5_array(df["CTX5"])

    for i in range(5):
        x[f"CTX5_{i}"] = ctx[:, i]

    x = x[FEATURE_NAMES]

    y = df["LABEL"].astype(np.int32).values

    meta = df[
        ["SAMPLE", "CHROM", "POS", "REF", "ALT", "VAR_TYPE", "LABEL"]
    ].copy()

    return x, y, meta


def safe_aucroc(y_true, pred):
    labels = np.unique(y_true)

    if len(labels) < 2:
        return None

    return float(roc_auc_score(y_true, pred))


def confusion_at_threshold(y_true, pred, threshold):
    yhat = (pred >= threshold).astype(np.int8)

    tn, fp, fn, tp = confusion_matrix(y_true, yhat, labels=[0, 1]).ravel()

    precision = precision_score(y_true, yhat, zero_division=0)
    recall = recall_score(y_true, yhat, zero_division=0)
    f1 = f1_score(y_true, yhat, zero_division=0)

    return {
        "threshold": float(threshold),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "predicted_positive": int(tp + fp),
    }


def make_threshold_table(y_true, pred):
    fixed_thresholds = [
        0.01,
        0.02,
        0.05,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50,
        0.60,
        0.70,
        0.80,
        0.90,
        0.95,
        0.975,
        0.99,
        0.995,
        0.999,
    ]

    rows = []

    for threshold in fixed_thresholds:
        row = confusion_at_threshold(y_true, pred, threshold)
        row["selection_rule"] = "fixed"
        rows.append(row)

    precision, recall, thresholds = precision_recall_curve(y_true, pred)

    if len(thresholds) > 0:
        f1_values = []

        for p, r in zip(precision[:-1], recall[:-1]):
            f1_values.append((2 * p * r / (p + r)) if (p + r) > 0 else 0.0)

        best_idx = int(np.argmax(f1_values))
        best_threshold = float(thresholds[best_idx])

        row = confusion_at_threshold(y_true, pred, best_threshold)
        row["selection_rule"] = "best_valid_f1"
        rows.append(row)

    for target_precision in [0.25, 0.50, 0.75, 0.90]:
        candidates = []

        for i, threshold in enumerate(thresholds):
            p = precision[i]
            r = recall[i]

            if p >= target_precision:
                candidates.append((r, threshold, p))

        if candidates:
            candidates.sort(reverse=True, key=lambda item: item[0])
            _, threshold, _ = candidates[0]

            row = confusion_at_threshold(y_true, pred, float(threshold))
            row["selection_rule"] = (
                f"max_recall_at_precision_ge_{target_precision}"
            )
            rows.append(row)

    table = pd.DataFrame(rows)

    table["threshold_rounded"] = table["threshold"].round(8)
    table = table.drop_duplicates(
        subset=["threshold_rounded", "selection_rule"]
    ).drop(columns=["threshold_rounded"])

    return table


def save_predictions(meta, pred, out_path):
    out = meta.copy()
    out["PRED_PROBA"] = pred.astype(np.float32)
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")


def save_feature_importance(booster, outdir, model_name):
    rows = []

    for importance_type in ["weight", "gain", "cover"]:
        scores = booster.get_score(importance_type=importance_type)

        for feature, value in scores.items():
            rows.append(
                {
                    "feature": feature,
                    "importance_type": importance_type,
                    "value": value,
                }
            )

    if rows:
        df = pd.DataFrame(rows)
        path = Path(outdir) / f"{model_name}.feature_importance.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"Saved feature importance: {path}")


def write_run_manifest(args, train_files, valid_files, test_files):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "train_files": train_files,
        "valid_files": valid_files,
        "test_files": test_files,
        "outdir": str(outdir),
        "seed": args.seed,
        "nthread": args.nthread,
        "num_boost_round": args.num_boost_round,
        "early_stopping_rounds": args.early_stopping_rounds,
        "models": args.models,
        "xgboost_version": xgb.__version__,
        "pandas_version": pd.__version__,
        "python": os.sys.version,
    }

    path = outdir / "training_manifest.json"

    with path.open("w") as handle:
        json.dump(manifest, handle, indent=2)

    print(f"Saved training manifest: {path}")


def train_and_eval(args, train_files, valid_files, test_files, want_types, model_name):
    print("\n" + "=" * 80)
    print(f"Training {model_name}")
    print(f"Variant types: {sorted(want_types)}")
    print("=" * 80)

    train_df = load_group(train_files, want_types)
    valid_df = load_group(valid_files, want_types)
    test_df = load_group(test_files, want_types)

    train_df = train_df.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)

    train_pos = int((train_df["LABEL"] == 1).sum())
    valid_pos = int((valid_df["LABEL"] == 1).sum())
    test_pos = int((test_df["LABEL"] == 1).sum())

    print("Train rows:", len(train_df), "positives:", train_pos)
    print("Valid rows:", len(valid_df), "positives:", valid_pos)
    print("Test rows :", len(test_df), "positives:", test_pos)

    if train_pos == 0:
        raise RuntimeError(f"No positives in training data for {model_name}")
    if valid_pos == 0:
        raise RuntimeError(f"No positives in validation data for {model_name}")
    if test_pos == 0:
        raise RuntimeError(f"No positives in test data for {model_name}")

    x_train, y_train, _ = make_matrix(train_df)
    x_valid, y_valid, meta_valid = make_matrix(valid_df)
    x_test, y_test, meta_test = make_matrix(test_df)

    del train_df, valid_df, test_df
    gc.collect()

    dtrain = xgb.DMatrix(x_train, label=y_train, feature_names=FEATURE_NAMES)
    dvalid = xgb.DMatrix(x_valid, label=y_valid, feature_names=FEATURE_NAMES)
    dtest = xgb.DMatrix(x_test, label=y_test, feature_names=FEATURE_NAMES)

    pos = int((y_train == 1).sum())
    neg = int((y_train == 0).sum())
    scale_pos_weight = neg / pos if pos > 0 else 1.0

    params = {
        "objective": "binary:logistic",
        "eval_metric": "aucpr",
        "max_depth": args.max_depth,
        "eta": args.eta,
        "subsample": args.subsample,
        "colsample_bytree": args.colsample_bytree,
        "min_child_weight": args.min_child_weight,
        "lambda": args.reg_lambda,
        "alpha": args.reg_alpha,
        "tree_method": args.tree_method,
        "nthread": args.nthread,
        "seed": args.seed,
        "scale_pos_weight": scale_pos_weight,
    }

    print("XGBoost parameters:")
    print(json.dumps(params, indent=2))

    booster = xgb.train(
        params=params,
        dtrain=dtrain,
        num_boost_round=args.num_boost_round,
        evals=[(dtrain, "train"), (dvalid, "valid")],
        early_stopping_rounds=args.early_stopping_rounds,
        verbose_eval=args.verbose_eval,
    )

    best_ntree_limit = booster.best_iteration + 1

    valid_pred = booster.predict(dvalid, iteration_range=(0, best_ntree_limit))
    test_pred = booster.predict(dtest, iteration_range=(0, best_ntree_limit))

    valid_aucpr = float(average_precision_score(y_valid, valid_pred))
    test_aucpr = float(average_precision_score(y_test, test_pred))

    valid_aucroc = safe_aucroc(y_valid, valid_pred)
    test_aucroc = safe_aucroc(y_test, test_pred)

    outdir = Path(args.outdir)

    valid_threshold_table = make_threshold_table(y_valid, valid_pred)
    valid_threshold_path = outdir / f"{model_name}.valid_thresholds.tsv"
    valid_threshold_table.to_csv(valid_threshold_path, sep="\t", index=False)

    selected_thresholds = valid_threshold_table[
        valid_threshold_table["selection_rule"] != "fixed"
    ][["selection_rule", "threshold"]].copy()

    fixed_selected = pd.DataFrame(
        {
            "selection_rule": ["fixed_0.50", "fixed_0.90", "fixed_0.95", "fixed_0.99"],
            "threshold": [0.50, 0.90, 0.95, 0.99],
        }
    )

    selected_thresholds = pd.concat(
        [selected_thresholds, fixed_selected], ignore_index=True
    )
    selected_thresholds = selected_thresholds.drop_duplicates(
        subset=["selection_rule", "threshold"]
    )

    test_rows = []

    for row in selected_thresholds.itertuples(index=False):
        metrics = confusion_at_threshold(y_test, test_pred, float(row.threshold))
        metrics["selection_rule"] = row.selection_rule
        test_rows.append(metrics)

    test_threshold_table = pd.DataFrame(test_rows)
    test_threshold_path = outdir / f"{model_name}.test_at_valid_selected_thresholds.tsv"
    test_threshold_table.to_csv(test_threshold_path, sep="\t", index=False)

    model_path = outdir / f"{model_name}.json"
    metrics_path = outdir / f"{model_name}.metrics.json"

    valid_names = "_".join(valid_files.keys())
    test_names = "_".join(test_files.keys())

    valid_pred_path = outdir / f"{model_name}.{valid_names}.valid_predictions.tsv.gz"
    test_pred_path = outdir / f"{model_name}.{test_names}.test_predictions.tsv.gz"

    booster.save_model(model_path)
    save_predictions(meta_valid, valid_pred, valid_pred_path)
    save_predictions(meta_test, test_pred, test_pred_path)
    save_feature_importance(booster, outdir, model_name)

    metrics = {
        "model_name": model_name,
        "variant_types": sorted(list(want_types)),
        "seed": args.seed,
        "best_iteration": int(booster.best_iteration),
        "best_score": float(booster.best_score),
        "nthread": args.nthread,
        "train_samples": list(train_files.keys()),
        "valid_samples": list(valid_files.keys()),
        "test_samples": list(test_files.keys()),
        "train_rows": int(len(y_train)),
        "train_pos": int((y_train == 1).sum()),
        "train_neg": int((y_train == 0).sum()),
        "valid_rows": int(len(y_valid)),
        "valid_pos": int((y_valid == 1).sum()),
        "valid_neg": int((y_valid == 0).sum()),
        "test_rows": int(len(y_test)),
        "test_pos": int((y_test == 1).sum()),
        "test_neg": int((y_test == 0).sum()),
        "scale_pos_weight": float(scale_pos_weight),
        "valid_aucpr": valid_aucpr,
        "test_aucpr": test_aucpr,
        "valid_aucroc": valid_aucroc,
        "test_aucroc": test_aucroc,
        "model_path": str(model_path),
        "valid_thresholds_path": str(valid_threshold_path),
        "test_thresholds_path": str(test_threshold_path),
        "valid_predictions_path": str(valid_pred_path),
        "test_predictions_path": str(test_pred_path),
        "features": FEATURE_NAMES,
        "params": params,
    }

    with metrics_path.open("w") as handle:
        json.dump(metrics, handle, indent=2)

    print("\nSaved model:", model_path)
    print("Saved metrics:", metrics_path)
    print("Saved validation predictions:", valid_pred_path)
    print("Saved test predictions:", test_pred_path)
    print("Saved validation threshold table:", valid_threshold_path)
    print("Saved test threshold table:", test_threshold_path)

    print("\nMetrics:")
    print(json.dumps(metrics, indent=2))

    del x_train, x_valid, x_test
    del y_train, y_valid, y_test
    del dtrain, dvalid, dtest
    del meta_valid, meta_test
    del valid_pred, test_pred
    gc.collect()


def main():
    args = parse_args()

    if args.nthread is None:
        args.nthread = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    train_files = sample_file_list_to_dict(args.train)
    valid_files = sample_file_list_to_dict(args.valid)
    test_files = sample_file_list_to_dict(args.test)

    check_input_files(train_files, "train")
    check_input_files(valid_files, "valid")
    check_input_files(test_files, "test")

    print("Output directory:", outdir)
    print("Train samples:", train_files)
    print("Valid samples:", valid_files)
    print("Test samples:", test_files)
    print("Python:", os.sys.version)
    print("xgboost:", xgb.__version__)
    print("pandas:", pd.__version__)

    write_run_manifest(args, train_files, valid_files, test_files)

    if args.models in ("both", "snp"):
        train_and_eval(
            args=args,
            train_files=train_files,
            valid_files=valid_files,
            test_files=test_files,
            want_types={"SNP"},
            model_name="xgb_snp_rescue_classifier",
        )

    if args.models in ("both", "indel"):
        train_and_eval(
            args=args,
            train_files=train_files,
            valid_files=valid_files,
            test_files=test_files,
            want_types={"INS", "DEL"},
            model_name="xgb_indel_rescue_classifier",
        )

    print("\nAll done.")


if __name__ == "__main__":
    main()
