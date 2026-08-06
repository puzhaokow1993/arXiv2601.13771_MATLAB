function [dl, bt] = circleg(R)
% Returns decomposed geometry of a circle with radius R

% Circle: [1; center_x; center_y; radius]
gd = [1; 0; 0; R];

% Names matrix for the geometry (must be char array)
ns = char('C1');
ns = ns';

% Set formula
sf = 'C1';

[dl, bt] = decsg(gd, sf, ns);
end
