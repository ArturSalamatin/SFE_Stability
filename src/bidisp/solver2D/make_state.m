function state = make_state(IC)
%MAKE_STATE Summary of this function goes here
%   Detailed explanation goes here

state = [IC.p(:); IC.c(:); IC.y(:); IC.G(:); IC.x(:)];

end

