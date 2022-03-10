function tree = soma_subtract(soma_tree, tree)

    %get a single coordinates for the soma (in this case the mean)
    mx=mean(soma_tree.X);
    my=mean(soma_tree.Y);
    mz=mean(soma_tree.Z);
    %subtract soma coordinate from rest to get the soma at the
    %orgin of the axis (for overlaying individuals cells later)
    if ~isempty(tree)
        tree.X=tree.X-mx;
        tree.Y=(tree.Y-my)*-1;
        tree.Z=tree.Z-mz;
    end
    
% 
%     com_tree.X=com_tree.X-mx;
%     com_tree.Y=(com_tree.Y-my)*-1;
%     com_tree.Z=com_tree.Z-mz;
% 
%     apical_tree.X=apical_tree.X-mx;
%     apical_tree.Y=(apical_tree.Y-my)*-1;
%     apical_tree.Z=apical_tree.Z-mz;
%     soma_tree.X=soma_tree.X-mx;
%     soma_tree.Y=soma_tree.Y-my;
%     soma_tree.Z=soma_tree.Z-mz;
%     if ~isempty(surface_tree)
%         surface_tree.X=surface_tree.X-mx;
%         surface_tree.Y=(surface_tree.Y-my)*-1;
%         surface_tree.Z=surface_tree.Z-mz;
%     end
%     if ~isempty(axon_tree)
%         axon_tree.X=axon_tree.X-mx;
%         axon_tree.Y=(axon_tree.Y-my)*-1;
%         axon_tree.Z=axon_tree.Z-mz;
%     end