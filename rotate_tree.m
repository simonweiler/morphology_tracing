function tree = rotate_tree(tree, rotation, y_offset)
    
    pts2d = [tree.X, tree.Y - y_offset]';
    rot_mat = [cos(rotation) -sin(rotation); sin(rotation) cos(rotation)];
    rot_pts = rot_mat * pts2d;
    tree.X = rot_pts(1, :)';
    tree.Y = rot_pts(2, :)' + y_offset;

end