%torsional defomration at centre of fixed-fixed beam calc - Roarks Stress and Strain table 10.3 (1g) 

T_0 = 100000*9.81*1; %Central torque
T_A = -T_0/2; %Fixed end reaction torque
T_B = T_A; %Fixed end reaction torque

t = RailLinkage.VLinks(1).tw; 
h = RailLinkage.VLinks(1).h{1}(0)-t;
b = RailLinkage.VLinks(1).w{1}(0)-t/2;

%Compute K (Table 10.1, 25 (and 24))
a_K = (RailLinkage.VLinks(1).h{1}(0)/2);
c_K = (RailLinkage.VLinks(1).w{1}(0))-t;
K_L1 = a_K*t^3*((1/3)-0.21*(t/a_K)*(1-t^4/(12*(a_K^4))));
K_L2 = c_K*t^3*((1/3)-0.105*(t/c_K)*(1-t^4/(192*(c_K^4))));
% alpha_K = 0.07;
% D_K = 0;
K_L = K_L1+K_L2;
K = 2*K_L;
G = RailLinkage.VLinks(1).G; %CONFIRM!
C_w = (h^2*b^3*t/12)*((2*h+3*b)/(h+6*b));
E =RailLinkage.VLinks(1).E;
Beta = sqrt(K*G/(C_w*E));
l = RailLinkage.VLinks(1).L;
theta_max = (T_0/(C_w*E*Beta^3))*((Beta*l/4)-tanh(Beta*l/4));

%% No torsional warping:
T_0 = 100000*9.81*1;

t = RailLinkage.VLinks(1).tw;

a_K = RailLinkage.VLinks(1).h{1}(0)/2;
c_K = RailLinkage.VLinks(1).w{1}(0)-t;

K_L1 = a_K*t^3*((1/3)-0.21*(t/a_K)*(1-t^4/(12*a_K^4)));
K_L2 = c_K*t^3*((1/3)-0.105*(t/c_K)*(1-t^4/(192*c_K^4)));

K_L = K_L1 + K_L2;
J = 2*K_L;         % Saint-Venant torsion constant

G = RailLinkage.VLinks(1).G;
L = RailLinkage.VLinks(1).L;

theta_max = T_0*L/(4*G*J);