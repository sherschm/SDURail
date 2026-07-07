function RigidLinks = makeRigidLink()
    RigidLinks = SorosimLink('empty'); %create an empty class element
    
    %To create a rigid link
    RigidLinks.jointtype='F';
    RigidLinks.linktype='r'; %rigid
    %RigidLinks.npie=2;
    
    CS = 'R'; %'C' for circular, 'R' for circular, 'E' for circular
    L = 10;
    Rho = 7850;

    Eta = 1e6;

    %C-section properties
    h  = @(X1)0.776; %height as function of X1
    w  = @(X1)1.140; %flange width
   
    %[M,cx] = CCrossSectProperties(L,Rho,h,w,tw,tf);
    [M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for rectangular cross section
    
    RigidLinks.ld{1}=L;
    RigidLinks.L= L;
    RigidLinks.CS= 'R';
    RigidLinks.r{1}= [];
    RigidLinks.h= {h};
    RigidLinks.w= {w};
    RigidLinks.a= [];
    RigidLinks.b= [];
    RigidLinks.cx= cx;
    RigidLinks.gi = {eye(4)};
    RigidLinks.gf = {eye(4)};
    RigidLinks.M= M;
    %RigidLinks.E = E ;
    RigidLinks.Rho= Rho;
    %RigidLinks.Poi = Poi;
    %RigidLinks.Eta = Eta;
    %RigidLinks.G=G;
    RigidLinks.Kj= [];
    RigidLinks.Dj= [];
    
    RigidLinks.n_l= 20;
    RigidLinks.n_r= 5; %should be 5 for rectangle, 9 for C-shaped X-sect
    RigidLinks.color= [0.9572 0.4854 0.8003];
    RigidLinks.alpha= 1;
    RigidLinks.CPF= false;
    RigidLinks.PlotFn= @(g)CustomShapePlot(g);
    RigidLinks.Lscale= 0.0947;
end

MoverCarriage = makeRigidLink();
%Now make full linkage object
MoverLinkage = SorosimLinkage(MoverCarriage);

