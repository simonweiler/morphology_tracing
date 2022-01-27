# morphology_tracing
code/scripts for reading out and analysing swc files of traced neurons

SW 160420 in TB lab updated on 220127 in TM lab
code to load and extract as well as save parameters morpholgy for cells that have basal and apical trees;
see as example https://www.biorxiv.org/content/10.1101/2021.10.13.464276v1.full 

dependencies:
1) TREES toolbox avialable at https://www.treestoolbox.org/ 
2) modified TREES scripts: cat_tree_sw, root_tree_sw, disscet_tree_sw, stats_tree_sw (on this repo)
3) uipickfiles get from https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids

%IMPORTANT: the dimensions should be ideally set/scaled correctly before tracing
%and exporting as swc files
