Note: All the python files only runs in the python 2.7 environment
The invokation of this script requires the following options:

Mandatory options:
--models: model path
--results_path: results path
--annot: annotation file path
--input_networks: input network file path
--prefix: prefix for result files
--goids: path to go id file
--protids: path to protein id file

Options with default values:
--valid_type: cv or temporal holdout, default is True(which means cv)
--epochs: iterations for trainning deep learning model, default 10
--batch_size: sample size for trainning each iteration, default 128
--alpha: propagation parameter, default 0.98
--select_arch: architecture for deepNF model, default[6*2500, 1200, 6*2500]

Optional options:
--topK: top K prediction score for protein go term pair
--cutoff_hi: cutting off go terms with more than higher bound of # of annotations
--cutoff_lo: cutting off go terms with less than lower bound of # of annotations
--prep: Option to preprocess the input network files

Depending on the prefix, if one model has already been made then new model will not be trainned.
If the input network file is already processed by this script, then remove the option --prep, it will automatically use the input network file as a standard input network file.  
