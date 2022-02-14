%% MODIFY ACCORDINGLY%
experiments_main_folder   = 'C:\';
% get directory for cell with swc files 
swc_list = uipickfiles('FilterSpec',experiments_main_folder);
%%%%
[morph_data] = morphology_readout(swc_list,true)