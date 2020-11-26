# CS4850
# Curriculum Based Chromosome Reconstruction: CBCR
## Van Hovenga 

### Usage: 
**Matlab**: To use, type in the terminal CBCR(input, curricula, alpha, gamma_1, gamma_2, learning_rate, max_iter, final_iter)
**Parameters**
- input: A string for the path of the input file.
- curricula: Integer. The number of curricula to be trained.
- alpha: Number between 0 and 1. The scaling factor for the trained data (values close to .5 are recommended). 
- gamma_1: Number between 0 and 1. The scaling factor for the first moment estimator in the adam optimizer (.9 is recommended).
- gamma_2: Number between 0 and 1. The scaling factor for the second moment estimator in the adam optimizer (.999 is recommended). 
- learning_rate: Learning rate for the optimizer (values less than .15 are recommended).
- max_iter: The maximumum iterations for training on the individual curricula.
- final_iter: The maximum iterationos for training on the whole data set after completion of curriculum training. 

### Input:
There are two possible input formats.
1. Tuple Input format(preferred) : A hi-C contact file, each line contains 3 numbers (separated by a space) of a contact, position_1 position_2 interaction_frequencies
2. Square Matrix Input format: The square matrix is a comma seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of equal-sized regions of a chromosome.

### Output: 
All outputs will be saved in a folder called *Scores* in the working directory. 
- output.log: A .log file that displays the optimal conversion factor for the trained structure along with the corresponding root mean-squared error, Pearson correlation distance, and Spearman correlation distance. 
- name_CONVERT_FACTOR=cfN=n.pdb: A .pdb file that shows the name of the imput data (name), the conversion factor for the corresponding structure (cf), and the number of curricula (n).  

