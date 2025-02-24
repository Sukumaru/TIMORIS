import sys
import os
import pytest
import pandas as pd
import numpy as np

# Ensure the 'src' directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from anndata import AnnData
from timoris import preprocess_data  # Import from updated __init__.py

@pytest.fixture
def synthetic_csv(tmp_path):
    """Create a synthetic CSV file with enough genes for testing."""
    file_path = tmp_path / "test_data.csv"
    df = pd.DataFrame(np.random.rand(50, 1000), columns=[f'Gene{i}' for i in range(1000)])  # Increased genes
    df.index.name = "cell"
    df.to_csv(file_path)
    return str(file_path)

@pytest.fixture
def synthetic_h5ad(tmp_path):
    """Create a synthetic .h5ad file with enough genes for testing."""
    file_path = tmp_path / "test_data.h5ad"
    adata = AnnData(np.random.rand(50, 1000))  # Increased genes
    adata.var_names = [f'Gene{i}' for i in range(1000)]
    adata.write_h5ad(file_path)
    return str(file_path)

def test_preprocess_h5ad(synthetic_h5ad):
    """Test preprocessing of .h5ad files."""
    adata = preprocess_data(synthetic_h5ad, file_type="h5ad", min_genes=10, min_cells=1)
    assert adata.n_obs > 0, "No cells left after preprocessing"
    assert adata.n_vars > 0, "No genes left after preprocessing"
    assert "highly_variable" in adata.var.columns

def test_preprocess_csv(synthetic_csv):
    """Test preprocessing of .csv files."""
    adata = preprocess_data(synthetic_csv, file_type="csv", min_genes=10, min_cells=1)
    assert adata.n_obs > 0, "No cells left after preprocessing"
    assert adata.n_vars > 0, "No genes left after preprocessing"
    assert "highly_variable" in adata.var.columns

if __name__ == "__main__":
    pytest.main()
