# Spectral matching algorithm
 Source code for the spectral matching algorithm, along with input data for the IN718 dataset.

# Required libraries
 The spectral matching algorithm requires numpy, scipy, and matplotlib (for test case plotting) installed in the local python environment.

# Running the test cases
 To run the test cases, call /path/to/local/python test.py (making sure the above libraries are installed). The three cases are
 1. Simple superposition of 3 cosines.
 2. Superposition of 50 cosines parametrized by random amplitude, period, and phase. Note that a seed is set so the output will be consistent.
 3. Width and depth data for IN718 melt pool variations.
 For each case, the input data is plotted in red, and the output data is plotted in royal blue. The PSD plots follow the same convention.

# Using the SpectralMatcher class
 To use the class for your own melt pool fluctuation generation, you need 
 1. Evenly spaced melt pool fluctuation data (say, from a high fidelity simulation, or optical microscopy of single tracks)
 2. The resolution of the fluctuation data, either in distance or time units (microns were used for development)
 In general, the worflow would be to
1. Import the spectral matching class from `spectral_matcher.py`
2. Instantiate the `SpectralMatcher` class with the experimental data sequence and resolution as inputs
3. Call the `gen_equivalent_fluctuation` method to get an output sequence
   Note that `return_all=True` can be used to get information on the spectral densities of the input and output to verify that they are equivalent. This can be seen in `test.py`.
