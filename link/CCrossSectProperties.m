function [M,cx] = CCrossSectProperties(L,Rho,h_in,b_in,tw,tf)
 % h  = height as function of X1
 % b  = flange width
 % tw = web thickness
 % tf = flange thickness
% ---- Section properties (constant) ----

h=h_in(0.0);
b=b_in(0.0);

A_web = tw * ( - 2*tf);
A_fl  = b * tf;
A     = A_web + 2*A_fl;

% Mass
mass = Rho * A * L;

% Centroid along beam
cx = L/2;

% ---- Section inertias (about section centroid) ----

% Web
Iy_web = (tw*(h - 2*tf)^3)/12;
Iz_web = ((h - 2*tf)*tw^3)/12;

% Flanges
Iy_fl = (b*tf^3)/12;
Iz_fl = (tf*b^3)/12;

y_fl = (h/2 - tf/2);

Iy_sec = Iy_web + 2*(Iy_fl + A_fl*y_fl^2);
Iz_sec = Iz_web + 2*(Iz_fl);

% ---- Beam integrals ----

% ∫ A(x-cx)^2 dx = A * L^3 / 12
Iy = Rho * ( Iy_sec * L + A * L^3 / 12 );
Iz = Rho * ( Iz_sec * L + A * L^3 / 12 );

% Torsion ( for C-section): torsional constant J
Ix = Rho * ( (Iy_sec + Iz_sec) * L );

 %J = (1/3)*(2*b*tf^3 + (h-2*tf)*tw^3);
 %Ix = Rho * ( J * L );

% ---- Mass matrix ----

M = [Ix 0 0 0 0 0;
     0 Iy 0 0 0 0;
     0 0 Iz 0 0 0;
     0 0 0 mass 0 0;
     0 0 0 0 mass 0;
     0 0 0 0 0 mass];

end