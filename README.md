# Covariance_comparison

This code is capable of calculating the difference between two compressed covariance matrices. It follows the procedure described in [Tassia & Valerio 2020](https://arxiv.org/pdf/2107.04211.pdf).

The compressed matrices were compressed using [MOPED]( https://doi.org/10.1046/j.1365-8711.2000.03692.x).

How to use:

```python
import numpy as np
from GetDiff import FindDiff

# Our covariances are in a .txt file, so we can import them as:
CompCovariance = np.loadtxt('CompressedCovariance.txt', dtype=np.float64)
CompGaussianCovariance = np.loadtxt('CompressedGaussianCovariance.txt', dtype=np.float64)

blah = FindDiff(CompGaussianCovariance, CompCovariance, sample_size=500)
# Choose a small step size for testing, but if you want proper results, then
# make sure this is > 1000. It will take longer, so be patient.

# Get the difference with the associated error (68%).
blah.diff
# Returns dictionary where:
# blah.diff['diag'] is the difference for the diagonal elements.
# blah.diff['corr] is the difference for the correlation matrix.

blah.chain
# Returns dictionary with a 2D numpy array with the first column corresponding to
# the chi2 value, and the second column to the respective value of the difference.
# blah.chain['diag']  for the diagonal elements.
# blah.chain['corr]  for the correlation matrix.
```

The plots are WIP.

Feel free to use, but please cite [Tassia & Valerio 2020](https://arxiv.org/pdf/2107.04211.pdf).
If you have comments, questions, or feedback, please contact Tassia (tassia.aferreira@gmail.com)
