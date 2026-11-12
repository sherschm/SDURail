function M = defineSpatialInertia(mass, Ic, c)

% Skew-symmetric matrix of CoM position
cx = [ ...
    0    -c(3)   c(2);
    c(3)  0     -c(1);
    -c(2)  c(1)   0    ];

M = [ ...
    Ic + mass*cx*cx',   mass*cx;
    -mass*cx,            mass*eye(3) ];

end

