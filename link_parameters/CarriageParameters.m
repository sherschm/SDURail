function RailLink = makeCarriage()
    RailLink = SorosimLink('empty'); %create an empty class element
    
    %To create a rigid link
    RailLink.jointtype='F';
    RailLink.linktype='r';
    RailLink.npie=1;
    
    CS = 'R'; %'C' for circular, 'R' for circular, 'E' for circular
    L = 1;

    %L = 0.001;
    %L = 25.415;
    Rho = 7850;
    %E = 69e9;
    %Poi = 0.33;
    %G =  E/(2*(1+Poi));
    %Eta = 1e6;
    % h = @(X1)0.776;
    % w = @(X1)1.140;
    % r = @(X1)0.03; %function of normalzied X. X1 in [0 1]
    % 
    %[M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for circular cross section
    
    mass = 100000+5362.126;
    cx = L/2;
    M = zeros(6,6);
    M(1,1) = 0.01;
    M(2,2) = 0.01;
    M(3,3) = 0.01;
    M(4,4) = mass;
    M(5,5) = mass;
    M(6,6) = mass;
    %C-section properties
    h  = @(X1)0.776; %height as function of X1
    w  = @(X1)1.140; %flange width
    tw = 0.05; %web thickness
    tf = 0.05; %flange thickness
   
    %[M,cx] = CCrossSectProperties(L,Rho,h,w,tw,tf);
     %[M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for rectangular cross section
    
    RailLink.ld=[];
    RailLink.L= L;
    RailLink.CS= 'R';
    RailLink.r{1}= [];
    RailLink.h= h;
    RailLink.w= w;
    RailLink.a= [];
    RailLink.b= [];
    RailLink.cx= cx;
    RailLink.gi = [eye(4)];
    %RailLink.gi = [eye(3) [0.0;1;0.0];0 0 0 1];
    RailLink.gf = [eye(4)];
    RailLink.M= M;
    %RailLink.E = E ;
    RailLink.Rho= Rho;
    %RailLink.Poi = Poi;
    %RailLink.Eta = Eta;
    %RailLink.G=G;
    RailLink.Kj= [];
    RailLink.Dj= [];
    
    RailLink.n_l= 2;
    RailLink.n_r= 5; %should be 5 for rectangle, 9 for C-shaped X-sect
    RailLink.color= [0.9572 0.4854 0.8003];
    RailLink.alpha= 1;
    RailLink.CPF= false;
    RailLink.PlotFn= @(g)CustomShapePlot(g);
    RailLink.Lscale= 0.0947;
end

%Carriage = SorosimLink();

Carriage = makeCarriage();
%Now make full linkage object

RailCarriage = SorosimLinkage(Carriage);

