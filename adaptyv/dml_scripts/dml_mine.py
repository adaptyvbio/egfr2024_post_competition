# Implementation of double ML for partial correlation analysis
# Three main approaches discussed:
# 1. Using pingouin's partial_corr
# 2. Using pingouin's mediation_analysis
# 3. Custom double ML implementation

# Import required packages
import polars as pl
from loguru import logger
import numpy as np
from tqdm import tqdm
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import KFold
import seaborn as sb
sb.set()

def get_expression_encodings(expr_series: pl.Series) -> dict:
    """
    Returns a dictionary of different encodings for expression values using Polars

    Args:
        expr_series: pl.Series containing expression values

    Returns:
        dict: Dictionary containing different encodings:
            - binary: pl.Series of 1s (high) and 0s (not high)
            - categorical: pl.Series of lowercase categories
            - onehot: pl.DataFrame with one-hot encoded columns
    """
    # Lowercase all values
    expr_lower = expr_series.str.to_lowercase()

    # Binary encoding (high vs rest)
    binary = (expr_lower == 'high').cast(pl.Int8)

    # One-hot encoding using Polars
    # First get unique categories
    unique_categories = expr_lower.unique().sort()

    # Create one-hot columns for each category
    onehot_cols = [
        (expr_lower == category).cast(pl.Int8).alias(f"expr_{category}")
        for category in unique_categories
    ]
    onehot = pl.DataFrame(onehot_cols)

    # Return dict of encodings
    return {
        'binary': binary,
        'categorical': expr_lower,
        'onehot': onehot
    }

def run_dml_analysis(
    df: pl.DataFrame,
    model_scores: list = ['esm_pll', 'normalized_esm_pll', 'i_plddt', 'pae_interaction'],
    binding_col: str = 'kd',
    expression_col: str = 'expression',
    output_path: str = "dml_results.csv"
) -> pl.DataFrame:
    """
    Run DML analysis for all combinations of model scores and expression encodings

    Args:
        df: pl.DataFrame containing model scores and expression data
        model_scores: list of column names containing model scores to analyze
        binding_col: name of column containing binding data (kd values)
        expression_col: name of column containing expression data
        output_path: path to save results CSV

    Returns:
        pl.DataFrame containing analysis results
    """
    # Prepare binding values (neg log10 of kd)
    binding_values = (
        -np.log10(
            df[binding_col]
            .cast(pl.Float64)
            .fill_null(1.0)  # Fill nulls with 1.0
            .to_numpy()
        )
    )

    # Get expression encodings as Polars objects
    expr_encodings = get_expression_encodings(df[expression_col])

    # Prepare results storage
    results = []

    # Calculate total iterations for progress bar
    total_iterations = len(model_scores) * len(expr_encodings)

    # Run all combinations
    with tqdm(total=total_iterations) as pbar:
        for score_col in model_scores:
            score_values = df[score_col].to_numpy()

            for encoding_name, expr_encoded in expr_encodings.items():
                if encoding_name == 'onehot':
                    # Run DML for each one-hot column
                    for col in expr_encoded.columns:
                        effect, stderr = dml_partial_effect(
                            score_values,
                            expr_encoded[col].to_numpy(),
                            binding_values
                        )
                        results.append({
                            'model_score': score_col,
                            'expression_encoding': f'{encoding_name}_{col}',
                            'partial_effect': effect,
                            'std_err': stderr
                        })
                else:
                    # Run DML for binary/categorical
                    effect, stderr = dml_partial_effect(
                        score_values,
                        expr_encoded.to_numpy(),
                        binding_values
                    )
                    results.append({
                        'model_score': score_col,
                        'expression_encoding': encoding_name,
                        'partial_effect': effect,
                        'std_err': stderr
                    })
                pbar.update(1)

    # Create results DataFrame using Polars
    results_df = pl.DataFrame(results)
    results_df.write_csv(output_path)
    logger.info(f"Results saved to {output_path}")

    return results_df

def analyze_expression_binding_dml(parquet_path: str):
    """
    Load parquet file and run DML analysis
    """
    logger.info("Loading parquet file from {}", parquet_path)
    df = pl.read_parquet(parquet_path)

    # Print columns
    logger.info("Columns in dataframe: {}", df.columns)

    # Print unique values in selected column and their counts
    logger.info("Unique values in 'selected' column: {}",
               df.group_by('selected').agg(pl.count()).sort('selected'))

    # Filter to keep both "Adaptyv selection" and "Top 100" rows
    logger.info("Filtering to keep Adaptyv selection and Top 100 rows...")
    df = df.filter(~(pl.col("selected") == "No"))

    # Filter out specific ID prefix
    logger.info("Filtering out specific ID prefix...")
    df = df.filter(~pl.col("id").str.starts_with("8a28fcc4-b172-4373-b524-78b2c1"))

    # Assert we have samples left
    n_samples = df.shape[0]
    assert n_samples > 0, "No samples left after filtering for selected==True"
    logger.info("Number of samples after filtering: {}", n_samples)

    # Fill nulls in kd column first
    df = df.with_columns(
        pl.col('kd').cast(pl.Float64).fill_null(1.0).alias('kd')
    )

    # Check for nulls in key columns
    null_counts = df.null_count()
    key_columns = ['expression', 'kd', 'esm_pll', 'normalized_esm_pll', 'i_plddt', 'pae_interaction']

    for col in key_columns:
        nulls = null_counts.get_column(col)[0]
        if nulls > 0:
            logger.info(f"Rows with null values in {col}:")
            logger.info(df.filter(pl.col(col).is_null()).select(['id', 'round', 'sequence', 'selected', col]))
            raise ValueError(f"Found {nulls} null values in column {col} after filtering for selected rows")

    logger.info("Running DML analysis...")
    results = run_dml_analysis(df)

    return results
# Basic implementation of double ML partial effects
def dml_partial_effect(M_scores, g_values, f_values, n_splits=5):
    """
    Estimate partial effect of M on f controlling for g using double ML
    Args:
        M_scores: np.array of model predictions
        g_values: np.array of conditioning variable
        f_values: np.array of target variable
    Returns:
        partial_effect: float, estimated effect
        std_err: float, standard error
    """
    # Convert g_values to numeric if it's categorical
    if isinstance(g_values[0], (str, np.str_)):
        # Create a mapping of unique values to integers
        unique_vals = np.unique(g_values)
        val_to_int = {val: i for i, val in enumerate(unique_vals)}
        g_values = np.array([val_to_int[val] for val in g_values])

    # First stage: residualize M and f w.r.t g
    kf = KFold(n_splits=n_splits)
    M_resid = np.zeros_like(M_scores)
    f_resid = np.zeros_like(f_values)

    for train_idx, test_idx in kf.split(g_values):
        # Fit nuisance models on training fold
        g_M = GradientBoostingRegressor().fit(
            g_values[train_idx].reshape(-1,1),
            M_scores[train_idx]
        )
        g_f = GradientBoostingRegressor().fit(
            g_values[train_idx].reshape(-1,1),
            f_values[train_idx]
        )

        # Predict and residualize on test fold
        M_resid[test_idx] = M_scores[test_idx] - g_M.predict(g_values[test_idx].reshape(-1,1))
        f_resid[test_idx] = f_values[test_idx] - g_f.predict(g_values[test_idx].reshape(-1,1))

    # Second stage: regress residualized f on residualized M
    partial_effect = (M_resid * f_resid).mean() / np.var(M_resid)
    std_err = np.sqrt(np.var(f_resid) / (len(M_resid) * np.var(M_resid)))

    return partial_effect, std_err

# More sophisticated implementation that returns calibrated predictor
def dml_estimate(
    M_scores: np.ndarray,    # Model predictions
    g_values: np.ndarray,    # Conditioning variable
    f_values: np.ndarray,    # Target variable
    n_splits: int = 5
) -> tuple[float, float, callable]:
    """
    Double/debiased ML estimation of partial effects between M and (f,g).
    Uses cross-fitting to avoid overfitting in nuisance estimation.
    Returns both the effect size and a function for test-time prediction.

    Key differences from naive implementation:
    1. Proper handling of first-stage estimation uncertainty
    2. Returns calibrated predictor for test-time use
    3. Accounts for cross-fitting in variance estimation
    """
    kf = KFold(n_splits=n_splits)
    M_resid = np.zeros_like(M_scores)
    f_resid = np.zeros_like(f_values)
    g_models = []  # Store models for test-time prediction

    # First stage: Cross-fitted residualization
    for train_idx, test_idx in kf.split(g_values):
        # Estimate E[M|g] and E[f|g]
        g_M = GradientBoostingRegressor().fit(
            g_values[train_idx].reshape(-1,1),
            M_scores[train_idx]
        )
        g_f = GradientBoostingRegressor().fit(
            g_values[train_idx].reshape(-1,1),
            f_values[train_idx]
        )

        # Store for test time
        g_models.append((g_M, g_f))

        # Get orthogonalized residuals
        M_resid[test_idx] = M_scores[test_idx] - g_M.predict(
            g_values[test_idx].reshape(-1,1)
        )
        f_resid[test_idx] = f_values[test_idx] - g_f.predict(
            g_values[test_idx].reshape(-1,1)
        )

    # Second stage: Efficient score estimation
    theta = (M_resid * f_resid).mean() / np.var(M_resid)

    # Variance estimation accounting for cross-fitting
    resid_sq = (f_resid - theta * M_resid)**2
    std_err = np.sqrt(np.mean(resid_sq) / (n_splits * np.var(M_resid)))

    # Create test-time predictor using averaged models
    def predict(M_new: np.ndarray) -> np.ndarray:
        """
        Apply estimated effect at test time.
        Only requires M scores, no access to f or g needed.
        """
        predictions = []
        for g_M, g_f in g_models:
            # Approximate E[f|M] using stored models
            pred = theta * (M_new - g_M.predict(g_M.predict(M_new).reshape(-1,1)))
            predictions.append(pred)
        return np.mean(predictions, axis=0)

    return theta, std_err, predict

def visualize_dml_results(results_df: pl.DataFrame, output_path: str = "dml_results.pdf"):
    """
    Create visualization of DML results with error bars and interpretation
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Convert to pandas for easier plotting with seaborn
    pdf = results_df.to_pandas()
    
    # Set up the plot style
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create grouped bar plot with error bars
    sns.barplot(
        data=pdf,
        x='model_score',
        y='partial_effect',
        hue='expression_encoding',
        alpha=0.8,
        ax=ax
    )
    
    # Add error bars manually
    for i, model in enumerate(pdf['model_score'].unique()):
        model_data = pdf[pdf['model_score'] == model]
        x_pos = i
        
        for j, (_, row) in enumerate(model_data.iterrows()):
            # Calculate offset for grouped bars
            width = 0.8 / len(model_data)  # 0.8 is default width in seaborn
            offset = (j - len(model_data)/2 + 0.5) * width
            
            ax.errorbar(
                x=x_pos + offset,
                y=row['partial_effect'],
                yerr=row['std_err'],
                fmt='none',
                color='black',
                capsize=3
            )

    # Customize the plot
    plt.xticks(rotation=45, ha='right')
    plt.title('Partial Effects of Model Scores on Binding\nControlling for Expression', pad=20)
    plt.xlabel('Model Score Type')
    plt.ylabel('Partial Effect (with ±1 std err)')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_path)
    logger.info(f"Plot saved to {output_path}")

    # Print interpretation
    print("\nInterpretation of DML Results:")
    print("------------------------------")
    print("The partial effect measures the predictive association between model scores")
    print("and binding strength (negative log10 KD), after controlling for expression level.")
    print("\nInterpretation by encoding type:")
    print("- binary: Association when controlling for high vs non-high expression")
    print("- categorical: Association when controlling for all expression categories")
    print("- onehot: Separate associations when controlling for each expression category")
    print("\nPositive partial effect means:")
    print("→ Higher model scores predict stronger binding (lower KD)")
    print("  after accounting for expression level variation")
    print("\nStandard errors (error bars) indicate uncertainty:")
    print("- Smaller error bars = more precise estimates")
    print("- Association is statistically reliable if error bars don't cross zero")

    # Add summary statistics
    significant_effects = pdf[abs(pdf['partial_effect']) > 2 * pdf['std_err']]
    if len(significant_effects) > 0:
        print("\nStatistically reliable associations (|effect| > 2*stderr):")
        for _, row in significant_effects.iterrows():
            print(f"- {row['model_score']} with {row['expression_encoding']}:")
            print(f"  Association strength: {row['partial_effect']:.3f} ± {row['std_err']:.3f}")

# Example usage:
if __name__ == "__main__":
    logger.info("Starting DML analysis pipeline")
    results = analyze_expression_binding_dml("all_submissions.parquet")
    visualize_dml_results(results)
    logger.info("Analysis complete")

