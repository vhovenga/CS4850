# CS4850
# Curriculum Based Chromosome Reconstruction: CBCR
## Van Hovenga 

### Usage 
**Matlab**: To use, type in the terminal CBCR(input, curricula, alpha, gamma_1, gamma_2, learning_rate, max_iter, final_iter)
**Parameters**
- input: A string for the path of the input file.
- curricula: Integer. The number of curricula to be trained.
- alpha: Number between 0 and 1. The scaling factor for the trained data. Values close to .5 are recommended. 
- gamma_1: Number between 0 and 1. The scaling factor for the first moment estimator in the adam optimizer. .9 is recommended.
- gamma_2: Number between 0 and 1. The scaling factor for the second moment estimator in the adam optimizer. .999 is recommended. 
- learning_rate: Learning rate for the optimizer. Values less than .15 are recommended.
- max_iter: The maximumum iterations for training on the individual curricula.
- final_iter: The maximum iterationos for training on the whole data set after completion of curriculum training. 

