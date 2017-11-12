# Notes on the MATLAB codes 

We provide the implementation of single-profile and multi-profile alignment methods, with and without using Gaussian process (GP) prior. Applications of the 
methods to the proteomic and glycomic data sets are demonstrated in the provided MATLAB codes: 

* Single-profile alignment (without GP prior): `bioinfo_proteo_sp.m`, `bioinfo_glyco_sp.m`
* Single-profile alignment (with GP prior): `bioinfo_proteo_gpsp.m`, `bioinfo_glyco_gpsp.m`
* Multi-profile alignment: (without GP prior): `bioinfo_proteo_mp4.m`, `bioinfo_glyco_mp4.m`
* Multi-profile alignment: (with GP prior): `bioinfo_proteo_gpmp4.m`, `bioinfo_glyco_gpmp4.m`

Please download associated MATLAB data matrices: 

* Base peak chromatograms for single-profile alignmnet: 
[proteomics](http://omics.georgetown.edu/alignLCMS/bpc_proteomics.mat), 
[glycomics](http://omics.georgetown.edu/alignLCMS/bpc_glycomics.mat)
* Binned chromatograms for multi-profile alignment: 
[proteomics](http://omics.georgetown.edu/alignLCMS/data_matrix_proteomics.mat),
[glycomics](http://omics.georgetown.edu/alignLCMS/data_matrix_glycomics.mat)

We use the [SIMA](http://hci.iwr.uni-heidelberg.de/MIP/Software/sima.php) model to perform the peak matching step

* Peak lists as input to the SIMA model: 
[proteomics](http://omics.georgetown.edu/alignLCMS/proteo_sima.zip), 
[glycomics](http://omics.georgetown.edu/alignLCMS/glyco_sima.zip)

For performance evaluation of the peak matching results:

* Evaluation scripts: `eval_proteo.m`, `eval_glyco.m`
* Ground-truth data: 
[proteomics](http://omics.georgetown.edu/alignLCMS/ground_truth_proteomics.mat),
[glycomics](http://omics.georgetown.edu/alignLCMS/ground_truth_glycomics.mat)

## Dependency

In order to run the MATLAB codes, you need to have the following utilities: 

### Attached functions

* `coda()` to compute MCQ values
* `priorRatioLn()` to compute prior ratio when not using Gaussian process prior
* `maxind()` to find local maxima in a sequence

### Functions from elsewhere 

* `inv_posdef()`, `randnorm()`, `scale_rows()`, `ndsum()` from Tom Minka's Lightspeed toolbox, available [here](http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/)
* `randraw()` from File Exchange at MATLAB Central, available [here](http://www.mathworks.com/matlabcentral/fileexchange/7309)
* `bsplinebasis()` from Scott Gaffney's CCToolbox, available [here](http://www.ics.uci.edu/~sgaffney/software/CCT/)
* `apcluster()` by Frey Lab, available [here](http://www.psi.toronto.edu/index.php?q=affinity%20propagation)
* GPML toolbox (v3.1) by Carl Edward Rasmussen and Hannes Nickisch, available [here](http://www.gaussianprocess.org/gpml/code/matlab/doc/index.html)

