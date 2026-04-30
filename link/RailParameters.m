function RailLink = makeRailLink()
    RailLink = SorosimLink('empty'); %create an empty class element
    
    %To create a rigid link
    RailLink.jointtype='F';
    RailLink.linktype='s';
    RailLink.npie=2;
    
    CS = 'R'; %'C' for circular, 'R' for circular, 'E' for circular
    L = 5.415;
    %L = 25.415;
    Rho = 1000;
    E = 6.9e8;
    Poi = 0.33;
    G =  E/(2*(1+Poi));
    Eta = 1e6;
    %h = @(X1)0.05;
    %w = @(X1)0.005;
    %r = @(X1)0.03; %function of normalzied X. X1 in [0 1]
    
    %[M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for circular cross section
    
    %C-section properties
    h  = @(X1)0.776; %height as function of X1
    w  = @(X1)1.140; %flange width
    tw = 0.05; %web thickness
    tf = 0.05; %flange thickness
   
    [M,cx] = CCrossSectProperties(L,Rho,h,w,tw,tf);
     %[M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for rectangular cross section
    
    RailLink.ld{1}=L;
    RailLink.L= L;
    RailLink.CS= 'R';
    RailLink.r{1}= [];
    RailLink.h= {h};
    RailLink.w= {w};
    RailLink.a= [];
    RailLink.b= [];
    RailLink.cx= cx;
    RailLink.gi = {eye(4)};
    RailLink.gf = {eye(4)};
    RailLink.M= M;
    RailLink.E = E ;
    RailLink.Rho= Rho;
    RailLink.Poi = Poi;
    RailLink.Eta = Eta;
    RailLink.G=G;
    RailLink.Kj= [];
    RailLink.Dj= [];
    
    RailLink.n_l= 25;
    RailLink.n_r= 5; %should be 5 for rectangle, 9 for C-shaped X-sect
    RailLink.color= [0.9572 0.4854 0.8003];
    RailLink.alpha= 1;
    RailLink.CPF= false;
    RailLink.PlotFn= @(g)CustomShapePlot(g);
    RailLink.Lscale= 0.0947;
end

RailLink = makeRailLink();
%Now make full linkage object
RailLinkage = SorosimLinkage(RailLink);

