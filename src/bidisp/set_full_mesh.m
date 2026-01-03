function mesh = set_full_mesh(N, params, delta)

left_mesh = set_left_mesh(N, params, delta);
right_mesh = set_right_mesh(N, params, delta);

mesh.N = left_mesh.N + right_mesh.N;
mesh.left = left_mesh;
mesh.right = right_mesh;

mesh.xBar = [left_mesh.xBar, right_mesh.xBar];
mesh.xBarL = left_mesh.xBarL;
mesh.xBarR = right_mesh.xBarR;
mesh.x = [left_mesh.x, right_mesh.x];

mesh.xi = [left_mesh.xi, right_mesh.xi];
mesh.xiL = left_mesh.xiL;
mesh.xiR = right_mesh.xiR;

% ids of all segments
mesh.segmIds = ...
    [1:(left_mesh.N-1), ...
    left_mesh.N+(1:(right_mesh.N-1))]; 

mesh.jumpIds = left_mesh.N;

end