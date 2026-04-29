%Function to compute boundary points of a cross section
%Last modified by Anup Teejo Mathew 28.11.2024
function [y,z] = computeBoundaryCBeamYZ(Link,X,varargin) % X varies from 0 to 1, varargin is division number (only for soft link)

    n_r = Link.n_r;
    
    j = varargin{1};
    h_fn  = Link.h{j};
    w_fn  = Link.w{j};
    t  = 0.05; %Link.t{j};

    h = h_fn(X);
    w = w_fn(X);

    % Half dimensions
    h2 = h/2;
    w2 = w/2;

    % Define C-shape (open on right side)
    y = [ ...
        -h2, -h2,  h2,  h2,  h2-t,  h2-t, ...
        -h2+t, -h2+t, -h2 ...
    ];

    z = [ ...
        -w2,  w2,  w2,  -w2,  -w2, w2-t, ...
        w2-t, -w2, -w2 ...
    ];

end