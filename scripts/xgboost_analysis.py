#!/usr/bin/env python3
"""
XGBoost Correlation, Feature Importance, and Feature Selection Analysis Script.

Performs correlation modeling using XGBoost to analyze relationships between
target variable and feature columns in CSV and Parquet files. Supports various
analysis methods including feature importance, regression, correlation analysis,
SHAP explanations, and report-only feature selection (without dropping columns).

Supported file formats:
- CSV (.csv)
- Parquet (.parquet, .pq) - requires pyarrow: pip install pyarrow
"""

import argparse
import sys
import signal
import logging
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pandas.api.types as ptypes

# Exit codes
STATUS_OK = 0
STATUS_ERROR = 1
STATUS_INTERRUPTED = 130

# Analysis methods
METHOD_FEATURE_IMPORTANCE = 'feature_importance'
METHOD_REGRESSION = 'regression'
METHOD_CORRELATION = 'correlation'
METHOD_SHAP = 'shap'
METHODS = [METHOD_FEATURE_IMPORTANCE, METHOD_REGRESSION, METHOD_CORRELATION, METHOD_SHAP]

# Feature-selection report methods
FEATURE_SELECTION_METHOD_IMPORTANCE = 'importance'
FEATURE_SELECTION_METHOD_PERMUTATION = 'permutation'
FEATURE_SELECTION_METHOD_MUTUAL_INFO = 'mutual_info'
FEATURE_SELECTION_METHOD_RFECV = 'rfecv'
FEATURE_SELECTION_METHODS = [
    FEATURE_SELECTION_METHOD_IMPORTANCE,
    FEATURE_SELECTION_METHOD_PERMUTATION,
    FEATURE_SELECTION_METHOD_MUTUAL_INFO,
    FEATURE_SELECTION_METHOD_RFECV,
]


def setup_logging(verbose=False, debug=False):
    """Configure logging with color formatting."""
    class ColorFormatter(logging.Formatter):
        COLORS = {
            'DEBUG': '\033[36m',
            'INFO': '\033[32m',
            'WARNING': '\033[33m',
            'ERROR': '\033[31m',
            'SUCCESS': '\033[92m',
            'RESET': '\033[0m'
        }

        def format(self, record):
            if record.levelname == 'INFO' and hasattr(record, 'success') and record.success:
                record.levelname = 'SUCCESS'

            color = self.COLORS.get(record.levelname, self.COLORS['RESET'])
            reset = self.COLORS['RESET']
            time_str = datetime.now().strftime('%H:%M:%S')
            return f"{time_str} - {color}{record.levelname}{reset} : {record.getMessage()}"

    logger = logging.getLogger()
    logger.handlers.clear()

    handler = logging.StreamHandler()
    handler.setFormatter(ColorFormatter())
    logger.addHandler(handler)

    if debug:
        logger.setLevel(logging.DEBUG)
    elif verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    return logger


def log_success(message: str):
    """Log a success message with special formatting."""
    logger = logging.getLogger()
    record = logging.LogRecord(
        name=logger.name, level=logging.INFO, pathname='', lineno=0,
        msg=message, args=(), exc_info=None
    )
    record.success = True
    logger.handle(record)


def signal_handler(signum, frame):
    """Handle interrupt signals gracefully."""
    print("\nOperation cancelled by user.")
    sys.exit(STATUS_INTERRUPTED)


def get_parser():
    """Get appropriate argument parser with color support if available."""
    try:
        from rich_argparse import ArgumentParser
        return ArgumentParser
    except ImportError:
        try:
            from argparse_color_formatter import ColorHelpFormatter

            class ColorArgumentParser(argparse.ArgumentParser):
                def __init__(self, *args, **kwargs):
                    kwargs.setdefault('formatter_class', ColorHelpFormatter)
                    super().__init__(*args, **kwargs)
            return ColorArgumentParser
        except ImportError:
            return argparse.ArgumentParser


class XGBoostAnalyzer:
    """XGBoost-based correlation, importance, and feature-selection analyzer."""

    def __init__(self, use_gpu=False, verbose=False, tree_method=None):
        self.use_gpu = use_gpu
        self.verbose = verbose
        self.logger = logging.getLogger()
        self.model = None
        self.feature_names = None
        self.tree_method = tree_method
        self.device = None
        self.encoding_maps = {}
        self._configure_device()

    def _configure_device(self):
        """Configure XGBoost device (CPU or GPU with Metal support)."""
        if self.tree_method is not None and self.tree_method != 'auto':
            if self.tree_method == 'gpu_hist':
                self.device = 'cuda'
                log_success(f"Tree method: {self.tree_method} (GPU)")
            else:
                self.device = 'cpu'
                log_success(f"Tree method: {self.tree_method} (CPU)")
            return

        if self.use_gpu:
            try:
                import platform
                if platform.system() == 'Darwin':
                    self.tree_method = 'hist'
                    self.device = 'cpu'
                    self.logger.info("Configuring XGBoost with optimized hist method (multi-core CPU)")
                    log_success("Optimized CPU mode enabled (hist tree method)")
                else:
                    self.tree_method = 'gpu_hist'
                    self.device = 'cuda'
                    self.logger.info("Configuring XGBoost for CUDA GPU acceleration")
                    log_success("GPU acceleration enabled (CUDA)")
            except Exception as exc:
                self.logger.warning(f"GPU initialization failed, falling back to CPU: {exc}")
                self.tree_method = 'hist'
                self.device = 'cpu'
        else:
            self.tree_method = 'hist'
            self.device = 'cpu'
            self.logger.info("Using CPU for XGBoost computation")

    def _get_base_params(self):
        """Get base XGBoost parameters."""
        params = {
            'tree_method': self.tree_method,
            'device': self.device,
            'objective': 'reg:squarederror',
            'eval_metric': 'rmse',
            'seed': 42,
            'verbosity': 0,
        }
        self.logger.debug(f"XGBoost params: tree_method={self.tree_method}, device={self.device}")
        return params

    def _train_quick_model(self, X, y, feature_names, num_rounds=200, params_override=None):
        """Train a lightweight booster for auxiliary tasks (e.g., feature selection)."""
        import xgboost as xgb

        params = self._get_base_params()
        params.update({
            'max_depth': 6,
            'eta': 0.1,
            'subsample': 0.9,
            'colsample_bytree': 0.9,
        })
        if params_override:
            params.update(params_override)

        dtrain = xgb.DMatrix(X, label=y, feature_names=feature_names)
        booster = xgb.train(params, dtrain, num_boost_round=num_rounds, verbose_eval=False)
        return booster

    def _get_sklearn_regressor(self):
        """Create an sklearn-compatible XGBoost regressor."""
        import xgboost as xgb

        estimator = xgb.XGBRegressor(
            n_estimators=400,
            max_depth=6,
            learning_rate=0.1,
            subsample=0.9,
            colsample_bytree=0.9,
            objective='reg:squarederror',
            tree_method=self.tree_method,
            device=self.device,
            reg_alpha=0,
            reg_lambda=1,
            random_state=42,
            n_jobs=-1,
            eval_metric='rmse',
            verbosity=0,
        )
        return estimator

    def save_model_directory(self, output_dir, method_name, results_df=None, additional_data=None):
        """Save model, metadata, and results to a directory."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Creating model directory: {output_path}")

        if self.model is not None:
            model_path = output_path / 'model.json'
            self.model.save_model(str(model_path))
            log_success(f"Model saved to: {model_path}")

        metadata = {
            'method': method_name,
            'timestamp': datetime.now().isoformat(),
            'use_gpu': self.use_gpu,
            'tree_method': self.tree_method,
            'device': self.device,
            'feature_names': self.feature_names,
            'num_features': len(self.feature_names) if self.feature_names else 0,
        }
        if self.encoding_maps:
            metadata['categorical_encodings'] = self.encoding_maps
        if additional_data:
            metadata.update(additional_data)

        metadata_path = output_path / 'metadata.json'
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        log_success(f"Metadata saved to: {metadata_path}")

        if self.feature_names:
            features_path = output_path / 'features.txt'
            with open(features_path, 'w') as f:
                f.write('\n'.join(self.feature_names))
            self.logger.info(f"Feature names saved to: {features_path}")

        if results_df is not None:
            results_path = output_path / 'results.csv'
            results_df.to_csv(results_path, index=False)
            log_success(f"Results saved to: {results_path}")

        print("\n" + "=" * 70)
        print(f"MODEL DIRECTORY CREATED: {output_path}")
        print("=" * 70)
        print("Contents:")
        for file in sorted(output_path.iterdir()):
            size = file.stat().st_size
            size_str = f"{size:,} bytes" if size < 1024 else f"{size/1024:.1f} KB"
            print(f"  {file.name:<20} {size_str}")
        print("=" * 70 + "\n")

        return output_path

    def save_feature_selection_report(self, results_df, report_path, metadata=None):
        """Persist feature-selection rankings to disk."""
        if results_df is None:
            return
        report_path = Path(report_path)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(report_path, index=False)
        log_success(f"Feature selection report saved to: {report_path}")

        if metadata:
            meta_path = report_path.with_suffix('.json')
            with open(meta_path, 'w') as f:
                json.dump(metadata, f, indent=2)

    def load_and_prepare_data(self, file_path, target_column, ignore_columns=None):
        """
        Load CSV or Parquet file and prepare features and target without dropping feature columns.
        All categorical/object columns are label-encoded, datetime columns become timestamps,
        and constant columns are retained (per request).
        """
        self.logger.info(f"Loading data from: {file_path}")
        file_path_lower = str(file_path).lower()
        is_parquet = file_path_lower.endswith('.parquet') or file_path_lower.endswith('.pq')

        try:
            if is_parquet:
                try:
                    import pyarrow.parquet as pq
                    self.logger.info("Reading Parquet file...")
                    try:
                        df = pd.read_parquet(file_path)
                    except Exception as exc:
                        if "Multiple matches" in str(exc) or "duplicate" in str(exc).lower():
                            self.logger.warning("Parquet file has duplicate columns, repairing...")
                            parquet_file = pq.ParquetFile(file_path)
                            table = parquet_file.read()
                            columns = list(table.column_names)
                            seen = {}
                            new_columns = []
                            for col in columns:
                                if col in seen:
                                    seen[col] += 1
                                    new_columns.append(f"{col}_duplicate_{seen[col]}")
                                else:
                                    seen[col] = 0
                                    new_columns.append(col)
                            df = table.to_pandas()
                            df.columns = new_columns
                            duplicates = sum(v for v in seen.values())
                            self.logger.info(f"Repaired {duplicates} duplicate column(s)")
                        else:
                            raise
                except ImportError:
                    self.logger.error("Parquet file detected but pyarrow not installed")
                    self.logger.error("Install with: pip install pyarrow")
                    raise ImportError("pyarrow library required for Parquet files")
            else:
                self.logger.info("Reading CSV file...")
                df = pd.read_csv(file_path)

            log_success(f"Successfully loaded {len(df):,} rows, {len(df.columns)} columns")
        except Exception as exc:
            file_type = "Parquet" if is_parquet else "CSV"
            self.logger.error(f"Error reading {file_type} file: {exc}")
            raise

        if target_column not in df.columns:
            self.logger.error(f"Target column '{target_column}' not found in file")
            self.logger.debug(f"Available columns: {list(df.columns)}")
            raise ValueError(f"Target column '{target_column}' not found")

        ignore_columns = [col.strip() for col in ignore_columns] if ignore_columns else []
        ignore_set = set(ignore_columns + [target_column])

        feature_cols = [col for col in df.columns if col not in ignore_set]
        if not feature_cols:
            raise ValueError("No feature columns available after filtering")

        self.logger.info(f"Target column: {target_column}")
        if ignore_columns:
            self.logger.info(f"Ignoring columns: {ignore_columns}")
        self.logger.info(f"Using {len(feature_cols)} feature columns")

        X = df[feature_cols].copy()
        y = df[target_column].copy()

        # Convert target to numeric early
        if not ptypes.is_numeric_dtype(y):
            self.logger.info("Converting target column to numeric...")
            y = pd.to_numeric(y, errors='coerce')

        self.encoding_maps = {}

        # Normalize feature dtypes without dropping columns
        for col in X.columns:
            series = X[col]

            if ptypes.is_bool_dtype(series):
                X[col] = series.astype(np.int8)
                continue

            if ptypes.is_numeric_dtype(series):
                X[col] = pd.to_numeric(series, errors='coerce')
                continue

            if ptypes.is_datetime64_any_dtype(series):
                ts = pd.to_datetime(series, utc=True, errors='coerce')
                X[col] = ts.view('int64') / 1e9
                self.logger.info(f"Converted datetime column '{col}' to UNIX seconds")
                continue

            # Fallback: label-encode categoricals/objects
            codes, uniques = pd.factorize(series.astype(str), sort=True)
            X[col] = codes.astype(np.float32)
            self.encoding_maps[col] = [None if pd.isna(val) else str(val) for val in uniques]
            self.logger.info(f"Label-encoded column '{col}' with {len(uniques)} categories")

        # Remove rows with missing target or features
        missing_mask = X.isna().any(axis=1) | y.isna()
        missing_count = int(missing_mask.sum())
        if missing_count > 0:
            self.logger.warning(f"Found {missing_count} rows with missing values (feature or target); removing rows")
            X = X.loc[~missing_mask]
            y = y.loc[~missing_mask]

        # Impute any residual NaNs with column medians (rare, e.g., entire column)
        if X.isna().any().any():
            self.logger.warning("Residual NaNs detected after row filtering; applying median imputation")
            for col in X.columns:
                if X[col].isna().any():
                    median_val = X[col].median()
                    if np.isnan(median_val):
                        median_val = 0.0
                    X[col] = X[col].fillna(median_val)

        # Handle infinite values
        inf_mask_X = np.isinf(X.values).any(axis=1)
        if inf_mask_X.any():
            count = int(inf_mask_X.sum())
            self.logger.warning(f"Found {count} rows with infinite feature values; removing")
            X = X.loc[~inf_mask_X]
            y = y.loc[~inf_mask_X]

        inf_mask_y = np.isinf(y.values)
        if inf_mask_y.any():
            count = int(inf_mask_y.sum())
            self.logger.warning(f"Found {count} rows with infinite target values; removing")
            X = X.loc[~inf_mask_y]
            y = y.loc[~inf_mask_y]

        MAX_FLOAT32 = 3.4e38
        extreme_mask = (np.abs(X.values) > MAX_FLOAT32).any(axis=1)
        if extreme_mask.any():
            count = int(extreme_mask.sum())
            self.logger.warning(f"Found {count} rows with extreme feature values (>{MAX_FLOAT32}); removing")
            X = X.loc[~extreme_mask]
            y = y.loc[~extreme_mask]

        if len(X) == 0:
            raise ValueError("No samples remain after preprocessing")

        feature_names = X.columns.tolist()
        self.feature_names = feature_names

        X_values = X.values.astype(np.float32, copy=False)
        y_values = y.values.astype(np.float32, copy=False)

        log_success(f"Data prepared: {len(X_values)} samples, {len(feature_names)} features (no columns dropped)")
        return X_values, y_values, feature_names

    def run_feature_selection(self, X, y, feature_names, methods, top_k=None,
                              output_path=None, subsample=None, permutation_repeats=5,
                              rfecv_step=1, rfecv_min_features=5, rfecv_folds=5):
        """
        Generate feature-selection rankings without altering the feature matrix.
        """
        if not methods:
            self.logger.info("No feature selection methods requested; skipping.")
            return None

        ordered_methods = []
        for method in methods:
            if method == 'all':
                ordered_methods = FEATURE_SELECTION_METHODS.copy()
                break
            if method not in FEATURE_SELECTION_METHODS:
                self.logger.warning(f"Skipping unsupported feature selection method: {method}")
                continue
            if method not in ordered_methods:
                ordered_methods.append(method)

        if not ordered_methods:
            self.logger.info("No valid feature selection methods after filtering; skipping.")
            return None

        X_fs = X
        y_fs = y
        if subsample and subsample > 0 and len(X) > subsample:
            rng = np.random.default_rng(42)
            idx = rng.choice(len(X), subsample, replace=False)
            X_fs = X[idx]
            y_fs = y[idx]
            self.logger.info(f"Subsampled {subsample} rows (from {len(X)}) for feature selection computations")

        results_df = pd.DataFrame({'feature': feature_names})
        rank_columns = []

        if FEATURE_SELECTION_METHOD_IMPORTANCE in ordered_methods:
            self.logger.info("Running gain-based feature importance for selection report")
            booster = self._train_quick_model(X_fs, y_fs, feature_names, num_rounds=200)
            gain_importance = booster.get_score(importance_type='gain')
            weight_importance = booster.get_score(importance_type='weight')
            results_df['xgb_gain'] = results_df['feature'].map(lambda f: gain_importance.get(f, 0.0))
            results_df['xgb_weight'] = results_df['feature'].map(lambda f: weight_importance.get(f, 0.0))
            results_df['importance_rank'] = results_df['xgb_gain'].rank(method='dense', ascending=False)
            rank_columns.append('importance_rank')

        if FEATURE_SELECTION_METHOD_PERMUTATION in ordered_methods:
            self.logger.info("Running permutation importance (sklearn)")
            from sklearn.inspection import permutation_importance

            estimator = self._get_sklearn_regressor()
            estimator.fit(X_fs, y_fs)
            perm = permutation_importance(
                estimator,
                X_fs,
                y_fs,
                n_repeats=permutation_repeats,
                random_state=42,
                n_jobs=-1
            )
            results_df['permutation_importance'] = perm.importances_mean
            results_df['permutation_std'] = perm.importances_std
            results_df['permutation_rank'] = results_df['permutation_importance'].rank(method='dense', ascending=False)
            rank_columns.append('permutation_rank')

        if FEATURE_SELECTION_METHOD_MUTUAL_INFO in ordered_methods:
            self.logger.info("Computing mutual information scores")
            from sklearn.feature_selection import mutual_info_regression

            mi_scores = mutual_info_regression(X_fs, y_fs, random_state=42)
            results_df['mutual_info'] = mi_scores
            results_df['mutual_info_rank'] = results_df['mutual_info'].rank(method='dense', ascending=False)
            rank_columns.append('mutual_info_rank')

        if FEATURE_SELECTION_METHOD_RFECV in ordered_methods:
            self.logger.info("Running RFECV (recursive feature elimination with CV)")
            from sklearn.feature_selection import RFECV
            from sklearn.model_selection import KFold

            min_features = max(1, min(rfecv_min_features, len(feature_names)))
            step = max(1, rfecv_step)
            folds = max(2, rfecv_folds)

            estimator = self._get_sklearn_regressor()
            cv = KFold(n_splits=folds, shuffle=True, random_state=42)
            rfecv = RFECV(
                estimator=estimator,
                step=step,
                min_features_to_select=min_features,
                cv=cv,
                scoring='neg_mean_squared_error',
                n_jobs=-1,
            )
            rfecv.fit(X_fs, y_fs)
            results_df['rfecv_rank'] = rfecv.ranking_.astype(int)
            results_df['rfecv_selected'] = rfecv.support_.astype(bool)
            rank_columns.append('rfecv_rank')

        if rank_columns:
            results_df['composite_rank'] = results_df[rank_columns].mean(axis=1)
        else:
            results_df['composite_rank'] = np.nan

        results_df = results_df.sort_values('composite_rank').reset_index(drop=True)
        results_df['recommended'] = False

        if top_k and top_k > 0:
            max_rank = results_df['composite_rank'].nsmallest(min(top_k, len(results_df))).max()
            results_df.loc[results_df['composite_rank'] <= max_rank, 'recommended'] = True

        print("\n" + "=" * 90)
        print("FEATURE SELECTION REPORT (columns retained; rankings only)")
        print("=" * 90)
        print(f"Methods: {', '.join(ordered_methods)}")
        print(f"Samples used: {len(X_fs)} (subsample={'yes' if len(X_fs) != len(X) else 'no'})")
        top_display = min(top_k or 20, len(results_df))
        columns_to_show = ['feature', 'composite_rank']
        optional_cols = [
            'xgb_gain', 'permutation_importance', 'mutual_info',
            'rfecv_rank', 'rfecv_selected', 'recommended'
        ]
        for col in optional_cols:
            if col in results_df.columns:
                columns_to_show.append(col)

        print(results_df[columns_to_show].head(top_display).to_string(index=False))
        print("=" * 90 + "\n")

        if output_path:
            metadata = {
                'methods': ordered_methods,
                'top_k': top_k,
                'samples_used': len(X_fs),
                'timestamp': datetime.now().isoformat(),
            }
            self.save_feature_selection_report(results_df, output_path, metadata)

        return results_df

    # Existing analysis methods (feature importance, correlation, regression, SHAP)
    # remain unchanged except for docstring updates referencing retained features.

    def feature_importance_analysis(self, X, y, feature_names, output_dir=None):
        import xgboost as xgb

        self.logger.info("Starting feature importance analysis")
        params = self._get_base_params()
        params.update({
            'max_depth': 6,
            'eta': 0.3,
            'subsample': 0.8,
            'colsample_bytree': 0.8,
        })

        dtrain = xgb.DMatrix(X, label=y, feature_names=feature_names)
        num_rounds = 100

        if self.verbose:
            print(f"\nTraining for {num_rounds} rounds (no validation set for feature importance)...")

        self.model = xgb.train(
            params,
            dtrain,
            num_boost_round=num_rounds,
            verbose_eval=False
        )

        log_success("Model training completed")
        importance_dict = self.model.get_score(importance_type='gain')

        results = []
        for feature in feature_names:
            importance = importance_dict.get(feature, 0.0)
            results.append({'feature': feature, 'importance': importance})

        results_df = pd.DataFrame(results).sort_values('importance', ascending=False)
        total_importance = results_df['importance'].sum()
        results_df['importance_pct'] = (results_df['importance'] / total_importance * 100) if total_importance > 0 else 0.0

        print("\n" + "=" * 80)
        print("FEATURE IMPORTANCE ANALYSIS RESULTS")
        print("=" * 80)
        print(f"Method: XGBoost (tree_method={self.tree_method})")
        print(f"Total features analyzed (none dropped): {len(feature_names)}")
        print(f"Training rounds: {num_rounds}")
        print("\nTop 20 Most Important Features:")
        print("-" * 80)
        print(f"{'Rank':<6} {'Feature':<40} {'Importance':<15} {'% Total':<10}")
        print("-" * 80)

        for rank, row in enumerate(results_df.head(20).itertuples(index=False), start=1):
            print(f"{rank:<6} {row.feature:<40} {row.importance:<15.4f} {row.importance_pct:<10.2f}")

        print("=" * 80)

        if output_dir:
            self.save_model_directory(
                output_dir,
                METHOD_FEATURE_IMPORTANCE,
                results_df,
                additional_data={'num_rounds': num_rounds}
            )

        return results_df

    def correlation_analysis(self, X, y, feature_names, output_dir=None):
        self.logger.info("Starting correlation analysis")
        pearson_correlations = []
        for i, feature in enumerate(feature_names):
            try:
                corr = np.corrcoef(X[:, i], y)[0, 1]
                pearson_correlations.append(corr if not np.isnan(corr) else 0.0)
            except Exception as exc:
                self.logger.warning(f"Error calculating correlation for {feature}: {exc}")
                pearson_correlations.append(0.0)

        importance_df = self.feature_importance_analysis(X, y, feature_names)
        results = []
        for i, feature in enumerate(feature_names):
            importance_row = importance_df[importance_df['feature'] == feature]
            xgb_importance = importance_row['importance'].values[0] if not importance_row.empty else 0.0
            results.append({
                'feature': feature,
                'pearson_correlation': pearson_correlations[i],
                'abs_pearson_correlation': abs(pearson_correlations[i]),
                'xgboost_importance': xgb_importance,
            })

        results_df = pd.DataFrame(results).sort_values('abs_pearson_correlation', ascending=False)

        print("\n" + "=" * 90)
        print("CORRELATION ANALYSIS RESULTS")
        print("=" * 90)
        print(f"Total features analyzed: {len(feature_names)}")
        print("\nTop 20 Correlations:")
        print("-" * 90)
        print(f"{'Rank':<6} {'Feature':<35} {'Pearson r':<15} {'|r|':<10} {'XGB Imp.':<12}")
        print("-" * 90)

        for rank, row in enumerate(results_df.head(20).itertuples(index=False), start=1):
            print(f"{rank:<6} {row.feature:<35} {row.pearson_correlation:<15.4f} "
                  f"{row.abs_pearson_correlation:<10.4f} {row.xgboost_importance:<12.4f}")

        print("=" * 90)

        if output_dir:
            self.save_model_directory(
                output_dir,
                METHOD_CORRELATION,
                results_df
            )

        return results_df

    def regression_analysis(self, X, y, feature_names, output_dir=None, test_size=0.2,
                            max_depth=6, learning_rate=0.1, n_estimators=500,
                            min_child_weight=1, gamma=0, subsample=0.8,
                            colsample_bytree=0.8, reg_alpha=0, reg_lambda=1,
                            early_stopping_rounds=50):
        import xgboost as xgb
        from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
        from sklearn.model_selection import train_test_split

        self.logger.info("Starting regression analysis")

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=42
        )

        params = self._get_base_params()
        params.update({
            'max_depth': max_depth,
            'eta': learning_rate,
            'min_child_weight': min_child_weight,
            'gamma': gamma,
            'subsample': subsample,
            'colsample_bytree': colsample_bytree,
            'reg_alpha': reg_alpha,
            'reg_lambda': reg_lambda,
        })

        dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=feature_names)
        dtest = xgb.DMatrix(X_test, label=y_test, feature_names=feature_names)

        if self.verbose:
            print("\nTraining Progress (showing every 25 rounds):")
            print(f"{'Round':<8} {'Train-RMSE':<12} {'Test-RMSE':<12} {'Test-R²':<10}")
            print("-" * 50)

        evals_result = {}

        class R2DisplayCallback(xgb.callback.TrainingCallback):
            def __init__(self, outer, display_freq=25):
                self.outer = outer
                self.display_freq = display_freq

            def after_iteration(self, model, epoch, evals_log):
                if not self.outer.verbose:
                    return False
                if epoch % self.display_freq == 0:
                    train_pred = model.predict(dtrain)
                    test_pred = model.predict(dtest)
                    train_r2 = r2_score(y_train, train_pred)
                    test_r2 = r2_score(y_test, test_pred)
                    train_rmse = evals_log['train']['rmse'][-1] if 'train' in evals_log else 0
                    test_rmse = evals_log['test']['rmse'][-1] if 'test' in evals_log else 0
                    print(f"{epoch:<8} {train_rmse:<12.6f} {test_rmse:<12.6f} {test_r2:<10.6f}")
                return False

        callbacks = [R2DisplayCallback(self)] if self.verbose else []

        self.model = xgb.train(
            params,
            dtrain,
            num_boost_round=n_estimators,
            evals=[(dtrain, 'train'), (dtest, 'test')],
            early_stopping_rounds=early_stopping_rounds,
            evals_result=evals_result,
            verbose_eval=False,
            callbacks=callbacks
        )

        best_iteration = self.model.best_iteration
        log_success(f"Model training completed (best iteration: {best_iteration})")

        y_pred_train = self.model.predict(dtrain)
        y_pred_test = self.model.predict(dtest)

        train_rmse = np.sqrt(mean_squared_error(y_train, y_pred_train))
        train_mae = mean_absolute_error(y_train, y_pred_train)
        train_r2 = r2_score(y_train, y_pred_train)

        test_rmse = np.sqrt(mean_squared_error(y_test, y_pred_test))
        test_mae = mean_absolute_error(y_test, y_pred_test)
        test_r2 = r2_score(y_test, y_pred_test)

        print("\n" + "=" * 70)
        print("REGRESSION ANALYSIS RESULTS")
        print("=" * 70)
        print(f"Method: XGBoost Regression (tree_method={self.tree_method})")
        print(f"Training samples: {len(X_train)}")
        print(f"Test samples: {len(X_test)}")
        print(f"Features retained: {len(feature_names)}")
        print(f"Best iteration: {best_iteration}")
        print("\nTraining Set Performance:")
        print(f"  RMSE: {train_rmse:.6f}")
        print(f"  MAE:  {train_mae:.6f}")
        print(f"  R²:   {train_r2:.6f}")
        print("\nTest Set Performance:")
        print(f"  RMSE: {test_rmse:.6f}")
        print(f"  MAE:  {test_mae:.6f}")
        print(f"  R²:   {test_r2:.6f}")
        print("=" * 70)

        predictions_df = pd.DataFrame({
            'actual': y_test,
            'predicted': y_pred_test,
            'error': y_test - y_pred_test,
            'abs_error': np.abs(y_test - y_pred_test),
        })

        if output_dir:
            self.save_model_directory(
                output_dir,
                METHOD_REGRESSION,
                predictions_df,
                additional_data={
                    'train_samples': len(X_train),
                    'test_samples': len(X_test),
                    'test_size': test_size,
                    'best_iteration': best_iteration,
                    'train_rmse': train_rmse,
                    'train_mae': train_mae,
                    'train_r2': train_r2,
                    'test_rmse': test_rmse,
                    'test_mae': test_mae,
                    'test_r2': test_r2,
                    'hyperparameters': {
                        'max_depth': max_depth,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'min_child_weight': min_child_weight,
                        'gamma': gamma,
                        'subsample': subsample,
                        'colsample_bytree': colsample_bytree,
                        'reg_alpha': reg_alpha,
                        'reg_lambda': reg_lambda,
                        'early_stopping_rounds': early_stopping_rounds,
                    }
                }
            )

        return {
            'train_metrics': {'rmse': train_rmse, 'mae': train_mae, 'r2': train_r2},
            'test_metrics': {'rmse': test_rmse, 'mae': test_mae, 'r2': test_r2},
            'model': self.model,
            'predictions': predictions_df,
        }

    def shap_analysis(self, X, y, feature_names, output_dir=None, max_samples=1000):
        try:
            import shap
        except ImportError:
            self.logger.error("SHAP library not installed. Install with: pip install shap")
            raise ImportError("SHAP library required for SHAP analysis")

        import xgboost as xgb

        self.logger.info("Starting SHAP analysis")

        if self.model is None:
            params = self._get_base_params()
            params.update({'max_depth': 6, 'eta': 0.3})
            dtrain = xgb.DMatrix(X, label=y, feature_names=feature_names)
            self.logger.info("Training model for SHAP analysis...")
            self.model = xgb.train(params, dtrain, num_boost_round=100, verbose_eval=False)

        if len(X) > max_samples:
            self.logger.info(f"Subsampling {max_samples} samples for SHAP calculation")
            indices = np.random.choice(len(X), max_samples, replace=False)
            X_sample = X[indices]
        else:
            X_sample = X

        self.logger.info("Computing SHAP values...")
        explainer = shap.TreeExplainer(self.model)
        shap_values = explainer.shap_values(X_sample)

        mean_shap_values = np.abs(shap_values).mean(axis=0)
        results = [{'feature': feature, 'mean_abs_shap': mean_shap_values[i]} for i, feature in enumerate(feature_names)]
        results_df = pd.DataFrame(results).sort_values('mean_abs_shap', ascending=False)

        total_shap = results_df['mean_abs_shap'].sum()
        results_df['shap_pct'] = (results_df['mean_abs_shap'] / total_shap * 100) if total_shap > 0 else 0.0

        print("\n" + "=" * 80)
        print("SHAP ANALYSIS RESULTS")
        print("=" * 80)
        print(f"Samples analyzed: {len(X_sample)}")
        print(f"Total features: {len(feature_names)}")
        print("\nTop 20 Most Important Features (by SHAP values):")
        print("-" * 80)
        print(f"{'Rank':<6} {'Feature':<40} {'Mean |SHAP|':<15} {'% Total':<10}")
        print("-" * 80)

        for rank, row in enumerate(results_df.head(20).itertuples(index=False), start=1):
            print(f"{rank:<6} {row.feature:<40} {row.mean_abs_shap:<15.6f} {row.shap_pct:<10.2f}")

        print("=" * 80)

        if output_dir:
            self.save_model_directory(
                output_dir,
                METHOD_SHAP,
                results_df,
                additional_data={
                    'samples_analyzed': len(X_sample),
                    'max_samples': max_samples,
                }
            )

        return results_df, shap_values


def main():
    """Main entry point for XGBoost analysis script."""
    signal.signal(signal.SIGINT, signal_handler)

    ParserClass = get_parser()
    parser = ParserClass(
        description="XGBoost-based correlation, feature importance, SHAP, and feature-selection analysis for CSV/Parquet data",
        epilog="Examples:\n"
               "  CSV:     python %(prog)s data.csv -t target_column -m correlation -v --gpu\n"
               "  Parquet: python %(prog)s data.parquet -t target_column -m feature_importance"
    )

    parser.add_argument("input", help="Input CSV or Parquet file path (.csv, .parquet, .pq)")
    parser.add_argument("-t", "--target-column", required=True, help="Name of the target column")
    parser.add_argument(
        "-m", "--method",
        choices=METHODS,
        default=METHOD_FEATURE_IMPORTANCE,
        help=f"Analysis method (default: {METHOD_FEATURE_IMPORTANCE})"
    )
    parser.add_argument("-o", "--output", help="Output directory path for model and results (optional)")
    parser.add_argument("-I", "--ignore-columns", help="Comma-separated list of column names to ignore")
    parser.add_argument("--gpu", action="store_true", help="Enable GPU acceleration (Mac Metal optimized)")
    parser.add_argument(
        "--tree-method",
        choices=['auto', 'exact', 'approx', 'hist', 'gpu_hist'],
        default='auto',
        help="XGBoost tree construction algorithm (default: auto)"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--test-size", type=float, default=0.2, help="Test set size for regression (default: 0.2)")
    parser.add_argument("--max-shap-samples", type=int, default=1000, help="Maximum samples for SHAP analysis (default: 1000)")

    # Hyperparameters for regression
    parser.add_argument("--max-depth", type=int, default=6, help="Maximum tree depth (default: 6)")
    parser.add_argument("--learning-rate", type=float, default=0.1, help="Learning rate/eta (default: 0.1)")
    parser.add_argument("--n-estimators", type=int, default=500, help="Number of boosting rounds (default: 500)")
    parser.add_argument("--min-child-weight", type=int, default=1, help="Minimum sum of instance weight in child (default: 1)")
    parser.add_argument("--gamma", type=float, default=0, help="Minimum loss reduction for split (default: 0)")
    parser.add_argument("--subsample", type=float, default=0.8, help="Subsample ratio of training data (default: 0.8)")
    parser.add_argument("--colsample-bytree", type=float, default=0.8, help="Subsample ratio of columns per tree (default: 0.8)")
    parser.add_argument("--reg-alpha", type=float, default=0, help="L1 regularization term (default: 0)")
    parser.add_argument("--reg-lambda", type=float, default=1, help="L2 regularization term (default: 1)")
    parser.add_argument("--early-stopping-rounds", type=int, default=50, help="Early stopping rounds (default: 50)")

    # Feature-selection options (report-only; no automatic dropping)
    parser.add_argument(
        "--feature-selection",
        nargs='+',
        choices=FEATURE_SELECTION_METHODS + ['all'],
        help="Generate report-only feature selection rankings (no automatic dropping). "
             "Choose one or multiple methods, or 'all'."
    )
    parser.add_argument(
        "--feature-selection-top-k",
        type=int,
        help="Highlight top-k ranked features in the report (columns are still retained in modeling)."
    )
    parser.add_argument(
        "--feature-selection-output",
        help="Optional CSV path for the feature selection report. "
             "If omitted, defaults to <output>/feature_selection_results.csv when --output is provided."
    )
    parser.add_argument(
        "--feature-selection-subsample",
        type=int,
        default=20000,
        help="Maximum rows used during feature selection (default: 20,000). Use 0 to disable subsampling."
    )
    parser.add_argument(
        "--permutation-repeats",
        type=int,
        default=5,
        help="Number of shuffles for permutation importance (default: 5)"
    )
    parser.add_argument(
        "--rfecv-step",
        type=int,
        default=1,
        help="Number of features removed per RFECV step (default: 1)"
    )
    parser.add_argument(
        "--rfecv-min-features",
        type=int,
        default=5,
        help="Minimum number of features RFECV is allowed to keep (default: 5)"
    )
    parser.add_argument(
        "--rfecv-folds",
        type=int,
        default=5,
        help="Number of cross-validation folds for RFECV (default: 5)"
    )

    args = parser.parse_args()

    logger = setup_logging(args.verbose, args.debug)
    logger.info("Starting XGBoost analysis")
    logger.debug(f"Arguments: {vars(args)}")

    try:
        import xgboost as xgb
        version = getattr(xgb, '__version__', 'unknown')
        logger.info(f"XGBoost version: {version}")
        if not hasattr(xgb, 'DMatrix') or not hasattr(xgb, 'train'):
            logger.error("XGBoost installation is incomplete or corrupted")
            sys.exit(STATUS_ERROR)
        log_success("XGBoost successfully imported and verified")
    except ImportError as exc:
        logger.error(f"XGBoost not installed: {exc}")
        logger.error("Install with: pip install xgboost")
        sys.exit(STATUS_ERROR)

    if args.method == METHOD_REGRESSION or args.feature_selection:
        try:
            import sklearn  # noqa: F401
            logger.debug("scikit-learn available for regression/feature-selection tasks")
        except ImportError:
            logger.error("scikit-learn required for regression or feature selection. Install with: pip install scikit-learn")
            sys.exit(STATUS_ERROR)

    if args.method == METHOD_SHAP:
        try:
            import shap  # noqa: F401
        except ImportError:
            logger.error("SHAP required for SHAP analysis. Install with: pip install shap")
            sys.exit(STATUS_ERROR)

    ignore_columns = []
    if args.ignore_columns:
        ignore_columns = [col.strip() for col in args.ignore_columns.split(',')]
        logger.debug(f"Parsed ignore columns: {ignore_columns}")

    tree_method = None if args.tree_method == 'auto' else args.tree_method
    analyzer = XGBoostAnalyzer(use_gpu=args.gpu, verbose=args.verbose, tree_method=tree_method)

    try:
        X, y, feature_names = analyzer.load_and_prepare_data(
            args.input,
            args.target_column,
            ignore_columns
        )
    except Exception as exc:
        logger.error(f"Failed to load and prepare data: {exc}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(STATUS_ERROR)

    feature_selection_methods = args.feature_selection
    if feature_selection_methods:
        subsample = None if args.feature_selection_subsample == 0 else args.feature_selection_subsample
        fs_output_path = args.feature_selection_output
        if fs_output_path is None and args.output:
            fs_output_path = str(Path(args.output) / "feature_selection_results.csv")

        try:
            analyzer.run_feature_selection(
                X,
                y,
                feature_names,
                methods=feature_selection_methods,
                top_k=args.feature_selection_top_k,
                output_path=fs_output_path,
                subsample=subsample,
                permutation_repeats=args.permutation_repeats,
                rfecv_step=args.rfecv_step,
                rfecv_min_features=args.rfecv_min_features,
                rfecv_folds=args.rfecv_folds
            )
        except Exception as exc:
            logger.error(f"Feature selection report failed: {exc}")
            if args.debug:
                import traceback
                traceback.print_exc()
            sys.exit(STATUS_ERROR)

    try:
        if args.method == METHOD_FEATURE_IMPORTANCE:
            analyzer.feature_importance_analysis(X, y, feature_names, args.output)
        elif args.method == METHOD_CORRELATION:
            analyzer.correlation_analysis(X, y, feature_names, args.output)
        elif args.method == METHOD_REGRESSION:
            analyzer.regression_analysis(
                X, y, feature_names,
                output_dir=args.output,
                test_size=args.test_size,
                max_depth=args.max_depth,
                learning_rate=args.learning_rate,
                n_estimators=args.n_estimators,
                min_child_weight=args.min_child_weight,
                gamma=args.gamma,
                subsample=args.subsample,
                colsample_bytree=args.colsample_bytree,
                reg_alpha=args.reg_alpha,
                reg_lambda=args.reg_lambda,
                early_stopping_rounds=args.early_stopping_rounds
            )
        elif args.method == METHOD_SHAP:
            analyzer.shap_analysis(X, y, feature_names, args.output, args.max_shap_samples)

        log_success(f"{args.method.replace('_', ' ').title()} analysis completed successfully")
    except Exception as exc:
        logger.error(f"Analysis failed: {exc}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(STATUS_ERROR)

    sys.exit(STATUS_OK)


if __name__ == "__main__":
    main()