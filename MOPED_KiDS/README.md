# MOPED for KiDS-1000 covariance matrix and dataset.

These files are used for compressing or transforming the [KiDS-1000 cosmic shear covariance matrix and dataset](https://arxiv.org/pdf/2007.15633), which can be found [here](http://kids.strw.leidenuniv.nl/sciencedata.php).

The compression is done using [MOPED](https://doi.org/10.1046/j.1365-8711.2000.03692.x) and thetransformation follows the procedure described in [Tassia & Valerio 2020](https://arxiv.org/pdf/2107.04211.pdf).

 To use these files, please download and install [cosmosis](https://bitbucket.org/joezuntz/cosmosis/wiki/Home) as well as [KCAP](https://github.com/KiDS-WL/kcap). You can then adjust the paths in moped_params.ini accordingly.

 The fiducial parameter values in moped_values.ini are set to the bestfit values found by KiDS-1000.

The compression code takes the 270x270 covariance matrix and compresses it to a 13x13 covariance matrix.

How to use:

```python
python get_moped.py
```

To perform a compression, please make sure transform = 'comp' in get_moped.py.
