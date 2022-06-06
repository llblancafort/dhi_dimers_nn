# dhi_dimers_nn
Scripts for DHI dimers Neural Networks paper
The data set is provided in file dimer_summary.dat, a three-column file containing the dimer names following the nomenclature from Ref. 1, and the Grel and Eexc endpoints for each dimer. Python scripts are provided for fingerprint generation (input_layer_generator.py) and for NN training and prediction for each endpoint (NN_G_rel.py and NN_E_exc.py).
Fingerprint generation. The input_layer_generation.py script generates five files with the MVR, MVS, QBF, QBB and QBS input sets.
Training and prediction. Run the NN_G_rel.py or NN_E_exc.py script passing the name of the input set file as argument, eg `python NN_G_rel.py MVS.out`. The NN_G_rel.py script generates the following files:
- NN_G_rel_names.txt: Names of dimers in the validation and tests sets.
- NN_G_rel.txt: Final value of training, validation and test loss.
- NN_G_rel_DFT_energies.txt: DFT energies of the dimers in the validation and test sets, following the order of NN_G_rel_names.txt.
- NN_G_rel_predict.txt: Predicted dimer energies for the validation and test sets, following the order of NN_G_rel_names.txt.
- NN_G_rel.csv: History of training and validation losses for each epoch.
- NN_G_rel.png: Plot of NN_G_rel.csv contained in the `images/NN_G_rel` directory.
The NN_E_exc.py script generates analogous files with `G_rel` replaced by `E_exc`.
