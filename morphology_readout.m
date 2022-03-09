
function [morph_data] = morphology_readout(swc_list, plot_cell, save_path)
%SW 160420 in TB lab updated on 220127 in TM lab
%code to load and extract as well as save parameters morpholgy for cells that have basal and apical trees;
% see as example https://www.biorxiv.org/content/10.1101/2021.10.13.464276v1.full 
%dependencies:
%1) TREES toolbox avilable at https://www.treestoolbox.org/ 
%2)%modified TREES scripts: cat_tree_sw, root_tree_sw, disscet_tree_sw,
%stats_tree_sw
%3) uipickfiles get from https://uk.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids

%IMPORTANT: the dimensions should be ideally set correctly before tracing
%and exporting as swc files

if nargin < 3
    save_path = '';
end
%% Find indexes of apical, basal and soma and check whether naimg of user is correct
for i=1:size(swc_list,2);
if contains(swc_list{1,i},'apical','IgnoreCase',true)==1;
 branchtype(i)=0;
elseif contains(swc_list{1,i},'basal','IgnoreCase',true)==1;
 branchtype(i)=1;
elseif contains(swc_list{1,i},'soma','IgnoreCase',true)==1;
 branchtype(i)=2;
elseif contains(swc_list{1,i},'surface','IgnoreCase',true)==1;
 branchtype(i)=3;
elseif contains(swc_list{1,i},'axon','IgnoreCase',true)==1;
 branchtype(i)=4;
else 
  disp('file names need to contain either apical, basal or soma');
end
end
apical=find(branchtype==0);
basal=find(branchtype==1);
soma=find(branchtype==2);
surface=find(branchtype==3);
axon=find(branchtype==4);
%% Extract parameters for apical and basal using stats_tree from TREES toolbox
%basal: all basal trees are first concatenated-> one basal tree
com_tree=[];
if length(basal)>1
tree1=[];tree2=[];
[tree1] = load_tree (swc_list{1,basal(1)}, 'swc'); 
[tree2] = load_tree (swc_list{1,basal(2)}, 'swc'); 
com_tree=cat_tree(tree1,tree2);
if length(basal)>2
    for i=1:length(basal)-2
        temp_tree=[];
%         try        
        temp_tree=load_tree (swc_list{1,basal(i+2)}, 'swc');
        com_tree=cat_tree_sw(com_tree,temp_tree);      
%         catch exception
%         throw(exception) 
%         end
    end
end
else 
 com_tree=load_tree (swc_list{1,basal(1)}, 'swc'); 
end
com_tree.name='basal_combined';


surface_tree = [];
if ~isempty(surface)
    [surface_tree] = load_tree (swc_list{1,surface(1)}, 'swc');
end

axon_tree = [];
if ~isempty(axon)
    [axon_tree] = load_tree (swc_list{1,axon(1)}, 'swc');
end

%% Morphology soma subtraction
[soma_tree] = load_tree (swc_list{1,soma(1)}, 'swc');
%get a single coordinates for the soma (in this case the mean)
mx=mean(soma_tree.X);
my=mean(soma_tree.Y);
mz=mean(soma_tree.Z);
%subtract soma coordinate from rest to get the soma at the
%orgin of the axis (for overlaying individuals cells later)
com_tree.X=com_tree.X-mx;
com_tree.Y=(com_tree.Y-my)*-1;
com_tree.Z=com_tree.Z-mz;
apical_tree.X=apical_tree.X-mx;
apical_tree.Y=(apical_tree.Y-my)*-1;
apical_tree.Z=apical_tree.Z-mz;
soma_tree.X=soma_tree.X-mx;
soma_tree.Y=soma_tree.Y-my;
soma_tree.Z=soma_tree.Z-mz;
morphology_traces={apical_tree;com_tree;soma_tree};

if ~isempty(surface)
    %         surface_tree.X = medfilt1(surface_tree.X, 10);
    %         surface_tree.Y = medfilt1(surface_tree.Y, 10);
    
    surfPts = [surface_tree.X surface_tree.Y];
    somaPt = [morphData.soma_stats.mx morphData.soma_stats.my];
    [k, dist] = dsearchn(surfPts, somaPt);
    nClosePts = 500;
    closestPt = surfPts(k);
    surfPtsClose = surfPts(k - nClosePts:k+nClosePts);
    c = polyfit(surfPtsClose(:, 1), surfPtsClose(:, 2), 1);
    m = c(1);
    surfAngleRad = atan(m);
    surfAngleDeg = rad2deg(surfAngleRad);
end


%stats see https://www.treestoolbox.org/downloads/TREES_manual.pdf
%page 75
%apical
for i=1:length(apical)
apical_tree=[];
[apical_tree] = load_tree (swc_list{1,apical(i)}, 'swc');
% Removed options '-x' here because we want sholl.
apical_stats=stats_tree_sw(apical_tree,[],[],'-w');close all;
end
%basal
% Removed option '-x' here because we want sholl.
basal_stats=stats_tree_sw(com_tree,[],[],'-w');close all;
basal_stats.gstats.basaltrees=length(basal);






if plot_cell
    % %plot the cell itself 
    big=figure;mon_pos=get(0,'MonitorPositions');set(gcf,'color','w', 'menubar','figure', 'position',[mon_pos(1,3)-1200 2 500 500]); % [left, bottom, width, height]
    figure(big);HP=plot_tree(apical_tree,[],[],[],[],'-3l');hold on;plot_tree(com_tree,[],[],[],[],'-3l');plot_tree(soma_tree,[],[],[],[],'-3l');
end
%save data in cd folder as morph_data (traces + stats)
cd(swc_list{1,apical(1)}(1:end-(length(apical_tree.name)+4)));
morph_data.traces=morphology_traces;
morph_data.traces_surface=surface_tree;
morph_data.traces_axon=axon_tree;

morph_data.apical_stats=apical_stats;
morph_data.basal_stats=basal_stats;
morph_data.soma_stats.mx=mx;
morph_data.soma_stats.my=my;
morph_data.soma_stats.mz=mz;

if isempty(save_path)
    save('morph_data','morph_data');
else
    save(save_path,'morph_data');
end

   
end
   
   
   
    
