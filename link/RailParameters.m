function RailLink = makeRailLink(twebb,tflange)
    RailLink = SorosimLink('empty'); %create an empty class element
    
    %To create a rigid link
    RailLink.jointtype='F';
    RailLink.linktype='s';
    RailLink.npie=2;
    
    CS = 'R'; %'C' for circular, 'R' for circular, 'E' for circular
    L = 5.415;
    %L = 25.415;
    Rho = 2700;
    E = 6.9e10;
    Poi = 0.33;
    G =  E/(2*(1+Poi));
    Eta = 1e6;

    %[M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for circular cross section
    
    %C-section properties FROM LSP CAD
    % h  = @(X1)0.776; %height as function of X1
    % w  = @(X1)1.140; %flange width
    % twebb = 0.05; %web thickness
    % tflange = 0.05; %flange thickness
    
    %C-section properties smaller example
    h  = @(X1)0.08; %height as function of X1
    w  = @(X1)0.1; %flange width
   
    [M,cx] = CCrossSectProperties(L,Rho,h,w,twebb,tflange);
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
       
    RailLink.n_l= 20;
    RailLink.n_r= 5; %should be 5 for rectangle, 9 for C-shaped X-sect
    RailLink.color= [0.9572 0.4854 0.8003];
    RailLink.alpha= 1;
    RailLink.CPF= false;
    RailLink.PlotFn= @(g)CustomShapePlot(g);
    RailLink.Lscale= 0.0947;
end
twebb = 0.003;
tflange= 0.003;
RailLink = makeRailLink(twebb,tflange);
%Now make full linkage object
RailLinkage = SorosimLinkage(RailLink);

