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
    """Create a synthetic CSV file for testing."""
    file_path = tmp_path / "test_data.csv"
    df = pd.DataFrame(np.random.rand(50, 100), columns=[f'Gene{i}' for i in range(100)])
    df.index.name = "cell"
    df.to_csv(file_path)
    return str(file_path)

@pytest.fixture
def synthetic_h5ad(tmp_path):
    """Create a synthetic .h5ad file for testing."""
    file_path = tmp_path / "test_data.h5ad"
    adata = AnnData(np.random.rand(50, 100))
    adata.var_names = [f'Gene{i}' for i in range(100)]
    adata.write_h5ad(file_path)
    return str(file_path)

def test_preprocess_h5ad(synthetic_h5ad):
    """Test preprocessing of .h5ad files."""
    adata = preprocess_data(synthetic_h5ad, file_type="h5ad")
    assert adata.n_obs > 0
    assert adata.n_vars > 0
    assert "highly_variable" in adata.var.columns

def test_preprocess_csv(synthetic_csv):
    """Test preprocessing of .csv files."""
    adata = preprocess_data(synthetic_csv, file_type="csv")
    assert adata.n_obs > 0
    assert adata.n_vars > 0
    assert "highly_variable" in adata.var.columns

if __name__ == "__main__":
    pytest.main()
