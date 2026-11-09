%Function to compute boundary points of a cross section
%Last modified by Anup Teejo Mathew 28.11.2024
function [y,z] = computeBoundaryYZ(Link,X,varargin) % X varies from 0 to 1, varargin is division number (only for soft link)

n_r = Link.n_r;

if Link.CS=='C'

    if Link.linktype == 'r'
        r_fn  = Link.r;
    else
        j = varargin{1}; %division number
        r_fn  = Link.r{j};
    end
    r     = r_fn(X);
    theta = linspace(0,2*pi,n_r);
    y     = r*sin(theta);
    z     = r*cos(theta);

elseif Link.CS=='R'
    
    if Link.linktype == 'r'
        h_fn  = Link.h;
        w_fn  = Link.w;
    else
        j = varargin{1}; %division number
        h_fn  = Link.h{j};
        w_fn  = Link.w{j};
    end

    h     = h_fn(X);
    w     = w_fn(X);
    y     = [h/2 -h/2 -h/2 h/2 h/2];
    z     = [w/2 w/2 -w/2 -w/2 w/2];

elseif Link.CS=='E'
    
    if Link.linktype == 'r'
        a_fn  = Link.a;
        b_fn  = Link.b;
    else
        j = varargin{1}; %division number
        a_fn  = Link.a{j};
        b_fn  = Link.b{j};
    end

    a     = a_fn(X);
    b     = b_fn(X);
    theta = linspace(0,2*pi,n_r);
    y     = a*sin(theta);
    z     = b*cos(theta);
    
elseif Link.CS == 'SDU'
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

end