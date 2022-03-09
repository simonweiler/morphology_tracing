
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
%% Load trees
% Load trees now so we can rotate them before statistics.

% Load basal
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
        temp_tree=load_tree (swc_list{1,basal(i+2)}, 'swc');
        com_tree=cat_tree_sw(com_tree,temp_tree);      
    end
end
else 
 com_tree=load_tree (swc_list{1,basal(1)}, 'swc'); 
end
com_tree.name='basal_combined';

% Load apical
[apical_tree] = load_tree (swc_list{1,apical(1)}, 'swc');

% Load soma
[soma_tree] = load_tree (swc_list{1,soma(1)}, 'swc');

% Load surface if it exists
surface_tree = [];
if ~isempty(surface)
    [surface_tree] = load_tree (swc_list{1,surface(1)}, 'swc');
end

% Load axon if it exists
axon_tree = [];
if ~isempty(axon)
    [axon_tree] = load_tree (swc_list{1,axon(1)}, 'swc');
end


%% Morphology soma subtraction
%get a single coordinates for the soma (in this case the mean)
soma_mx=mean(soma_tree.X);
soma_my=mean(soma_tree.Y);
soma_mz=mean(soma_tree.Z);
% Save these for later, might be useful.
morph_data.soma_stats.mx=soma_mx;
morph_data.soma_stats.my=soma_my;
morph_data.soma_stats.mz=soma_mz;
%subtract soma coordinate from rest to get the soma at the
%orgin of the axis (for overlaying individuals cells later)
com_tree.X=com_tree.X-soma_mx;
com_tree.Y=(com_tree.Y-soma_my)*-1;
com_tree.Z=com_tree.Z-soma_mz;
apical_tree.X=apical_tree.X-soma_mx;
apical_tree.Y=(apical_tree.Y-soma_my)*-1;
apical_tree.Z=apical_tree.Z-soma_mz;
soma_tree.X=soma_tree.X-soma_mx;
%soma_tree.Y=soma_tree.Y-my; % why was this not flipped?
soma_tree.Y=(soma_tree.Y-soma_my)*-1;
soma_tree.Z=soma_tree.Z-soma_mz;
if ~isempty(surface_tree)
    surface_tree.X=surface_tree.X-soma_mx;
    surface_tree.Y=(surface_tree.Y-soma_my)*-1;
    surface_tree.Z=surface_tree.Z-soma_mz;
end
if ~isempty(axon_tree)
    axon_tree.X=axon_tree.X-soma_mx;
    axon_tree.Y=(axon_tree.Y-soma_my)*-1;
    axon_tree.Z=axon_tree.Z-soma_mz;
end
% mx=mean(soma_tree.X)
% my=mean(soma_tree.Y)


%% Rotate everything relative to the surface.
if ~isempty(surface_tree)
    
    surfPts = [surface_tree.X surface_tree.Y];
    nSurfPts = length(surface_tree.X);
    somaPt = [0 0]; % Assume soma has been centred to zero, see above.

    [k, surf_dist] = dsearchn(surfPts, somaPt);
    % todo this is hard coded pixels but could be a distance parameter if
    % the dpi is known.
    nClosePts = 500;
    closestPt = surfPts(k, :);
    indexPtStart = k - nClosePts;
    indexPtEnd = k + nClosePts;
    if indexPtStart < 1; indexPtStart = 1; end
    if indexPtEnd > nSurfPts; indexPtEnd = nSurfPts; end
    surfPtsClose = surfPts(indexPtStart:indexPtEnd, :);
    c = polyfit(surfPtsClose(:, 1), surfPtsClose(:, 2), 1);
    surfFitGrad = c(1);
    surfFitOffset = 0; %c(2);

    surfAngleRad = atan(-surfFitGrad);
    if closestPt(2) < 0
        % Surface is bellow soma, so needs to rotate the other way
        surfAngleRad = pi + surfAngleRad;
    end
    
    surfAngleDeg = rad2deg(surfAngleRad);

    apical_tree = rotate_tree(apical_tree, surfAngleRad, surfFitOffset);
    com_tree = rotate_tree(com_tree, surfAngleRad, surfFitOffset);
    soma_tree = rotate_tree(soma_tree, surfAngleRad, surfFitOffset);
    surface_tree = rotate_tree(surface_tree, surfAngleRad, surfFitOffset);
    if ~isempty(axon_tree)
        axon_tree = rotate_tree(axon_tree, surfAngleRad, surfFitOffset);
    end

end



%% Extract parameters for apical and basal using stats_tree from TREES toolbox
% Apical
% Removed option '-x' here because we want sholl.
apical_stats=stats_tree_sw(apical_tree,[],[],'-w');close all;
% Basal
% Removed option '-x' here because we want sholl.
basal_stats=stats_tree_sw(com_tree,[],[],'-w');close all;
basal_stats.gstats.basaltrees=length(basal);

if plot_cell
    % %plot the cell itself 
    big=figure;mon_pos=get(0,'MonitorPositions');set(gcf,'color','w', 'menubar','figure', 'position',[mon_pos(1,3)-1200 2 500 500]); % [left, bottom, width, height]
    figure(big);HP=plot_tree(apical_tree,[],[],[],[],'-3l');hold on;plot_tree(com_tree,[],[],[],[],'-3l');plot_tree(soma_tree,[],[],[],[],'-3l');
end

morphology_traces={apical_tree;com_tree;soma_tree};
morph_data.traces=morphology_traces;
morph_data.traces_opt.surface=surface_tree;
morph_data.traces_opt.axon=axon_tree;

morph_data.apical_stats=apical_stats;
morph_data.basal_stats=basal_stats;
if ~isempty(surface_tree)
    morph_data.surface_stats.dist_soma = surf_dist;
    morph_data.surface_stats.angle_soma_deg = surfAngleDeg;
else
    morph_data.surface_stats=NaN;
    morph_data.surface_stats.angle_soma_rad = Nan;
end

if isempty(save_path)
    %save data in cd folder as morph_data (traces + stats)
    cd(swc_list{1,apical(1)}(1:end-(length(apical_tree.name)+4)));
    save('morph_data','morph_data');
else
    save(save_path,'morph_data');
end

   
end
   
   
   
    
