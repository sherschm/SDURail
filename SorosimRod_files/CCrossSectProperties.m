% function [M,cx] = CCrossSectProperties(L,Rho,h_in,b_in,tw,tf)
% % h  = height as function of s
% % b  = flange width as function of s
% % tw = web thickness
% % tf = flange 5thickness
% % ---- Section properties (constant) ----
% 
% h=h_in(0.0);
% b=b_in(0.0);
% 
% A_web = tw * ( h - 2*tf);
% A_fl  = b * tf;
% A     = A_web + 2*A_fl;
% 
% % Mass
% mass = Rho * A * L;
% 
% % Centroid along beam
% cx = L/2;
% 
% % ---- Section inertias (about section centroid) ----
% 
% % Web
% Iy_web = (tw*(h - 2*tf)^3)/12;
% Iz_web = ((h - 2*tf)*tw^3)/12;
% 
% % Flanges
% Iy_fl = (b*tf^3)/12;
% Iz_fl = (tf*b^3)/12;
% 
% y_fl = (h/2 - tf/2);
% 
% Iy_sec = Iy_web + 2*(Iy_fl + A_fl*y_fl^2);
% Iz_sec = Iz_web + 2*(Iz_fl);
% 
% % ---- Beam integrals ----
% 
% % ∫ A(x-cx)^2 dx = A * L^3 / 12
% Iy = Rho * ( Iy_sec * L + A * L^3 / 12 );
% Iz = Rho * ( Iz_sec * L + A * L^3 / 12 );
% 
% % Torsion ( for C-section): torsional constant J
% Ix = Rho * ( (Iy_sec + Iz_sec) * L );
% 
% %J = (1/3)*(2*b*tf^3 + (h-2*tf)*tw^3);
% %Ix = Rho * ( J * L );
% 
% % ---- Mass matrix ----
% 
% M = [Ix 0 0 0 0 0;
%     0 Iy 0 0 0 0;
%     0 0 Iz 0 0 0;
%     0 0 0 mass 0 0;
%     0 0 0 0 mass 0;
%     0 0 0 0 0 mass];
% 
% end


function [M,cx,A,Iy_sec,Iz_sec,J] = ...
    CCrossSectProperties(L,Rho,h_in,b_in,tw,tf)

% C-section (channel section) rigid-body mass properties
%
% h  = overall depth
% b  = flange width
% tw = web thickness
% tf = flange thickness
%
% Returns:
%   M      = rigid-body mass matrix at COM
%   cx     = COM location along beam axis
%   A      = cross-sectional area
%   Iy_sec = centroidal second moment (about y-axis)
%   Iz_sec = centroidal second moment (about z-axis)
%   J      = Saint-Venant torsion constant

%% Section dimensions

h = h_in(0.0);
b = b_in(0.0);

hw = h - 2*tf;      % web height

%% Areas

A_web = tw*hw;
A_fl  = b*tf;

A = A_web + 2*A_fl;

%% Section centroid (z-direction)

z_web = tw/2;
z_fl = tw + (b - tw)/2;

zc = ( ...
    A_web*z_web ...
    + 2*A_fl*z_fl ) / A;

%% Centroidal Iy

Iy_web = tw*hw^3/12;

Iy_fl_local = b*tf^3/12;
y_fl = h/2 - tf/2;

Iy_sec = Iy_web ...
    + 2*(Iy_fl_local + A_fl*y_fl^2);

%% Centroidal Iz

Iz_web_local = hw*tw^3/12;
Iz_fl_local  = tf*b^3/12;

Iz_sec = ...
    (Iz_web_local + A_web*(z_web-zc)^2) ...
    + 2*(Iz_fl_local + A_fl*(z_fl-zc)^2);

%% Saint-Venant torsion constant
% Thin-walled open-section approximation

J = (1/3)*( ...
    2*b*tf^3 + ...
    hw*tw^3 );

%% Beam mass

mass = Rho*A*L;

%% Beam COM location

cx = L/2;

%% Rigid-body inertias about beam COM

Ixx = Rho*J*L;

Iyy = Rho*( ...
    Iy_sec*L + ...
    A*L^3/12 );

Izz = Rho*( ...
    Iz_sec*L + ...
    A*L^3/12 );

%% Rigid-body mass matrix

M = diag([ ...
    Ixx ...
    Iyy ...
    Izz ...
    mass ...
    mass ...
    mass ]);

end