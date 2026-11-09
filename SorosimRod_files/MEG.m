%This function computes Ms, Es, and Gs for a given Xs for the jth division
%of a link (15.12.2022)

function [Ms,Es,Gs]= MEG(Link,j,Xs)
np = length(Xs);
%j is division number

if Link.CS=='R'
    h_fn     = Link.h{j};
    w_fn     = Link.w{j};

    %updating:
    h_nGauss = zeros(np,1);
    w_nGauss = zeros(np,1);
    for ii=1:np
        h_nGauss(ii) = h_fn(Xs(ii));
        w_nGauss(ii) = w_fn(Xs(ii));
    end

    Iy_p = (1/12)*h_nGauss.*(w_nGauss.^3);
    Iz_p = (1/12)*(h_nGauss.^3).*w_nGauss;
    Ix_p = Iy_p+Iz_p;
    A_p  = h_nGauss.*w_nGauss;

elseif Link.CS=='C' %Circular
    r_fn     = Link.r{j};

    %updating:
    r_nGauss = zeros(np,1);
    for ii=1:np
        r_nGauss(ii) = r_fn(Xs(ii));
    end
    Iy_p = (pi/4)*r_nGauss.^4;
    Iz_p = Iy_p;
    Ix_p = Iy_p+Iz_p;
    A_p  = pi*r_nGauss.^2;

elseif Link.CS=='E'
    a_fn     = Link.a{j};
    b_fn     = Link.b{j};

    %updating:
    a_nGauss = zeros(np,1);
    b_nGauss = zeros(np,1);
    for ii=1:np
        a_nGauss(ii) = a_fn(Xs(ii));
        b_nGauss(ii) = b_fn(Xs(ii));
    end

    Iy_p = (pi/4)*a_nGauss.*(b_nGauss.^3);
    Iz_p = (pi/4)*(a_nGauss.^3).*b_nGauss;
    Ix_p = Iy_p+Iz_p;
    A_p  = pi*a_nGauss.*b_nGauss;

elseif Link.CS=='SDU'   % Channel ("C") section

    h_fn  = Link.h{j};
    w_fn  = Link.w{j};
    tw = Link.tw;
    tf = Link.tf;

    A_p  = zeros(np,1);
    Iy_p = zeros(np,1);
    Iz_p = zeros(np,1);
    Ix_p = zeros(np,1);

    for ii = 1:np

        h  = h_fn(Xs(ii));
        w  = w_fn(Xs(ii));

        hw = h - 2*tf;

        % Areas
        A_web = tw*hw;
        A_fl  = w*tf;

        A = A_web + 2*A_fl;

        % Centroid (z-direction)
        z_web = tw/2;
        z_fl  = tw + (w-tw)/2;

        zc = ( ...
            A_web*z_web ...
            + 2*A_fl*z_fl ) / A;

        % Centroidal Iy
        Iy_web = tw*hw^3/12;

        Iy_fl_local = w*tf^3/12;
        y_fl = h/2 - tf/2;

        Iy = Iy_web ...
            + 2*(Iy_fl_local + A_fl*y_fl^2);

        % Centroidal Iz
        Iz_web_local = hw*tw^3/12;
        Iz_fl_local  = tf*w^3/12;

        Iz = ...
            (Iz_web_local + A_web*(z_web-zc)^2) ...
            + 2*(Iz_fl_local + A_fl*(z_fl-zc)^2);

        % Saint-Venant torsion constant
        J = (1/3)*( ...
            2*w*tf^3 ...
            + hw*tw^3 );

        A_p(ii)  = A;
        Iy_p(ii) = Iy;
        Iz_p(ii) = Iz;
        Ix_p(ii) = J;
    end
    
end

Ms = zeros(6*np,6); %inertia
Es = zeros(6*np,6); %stiffness
Gs = zeros(6*np,6); %damping

Rho  = Link.Rho;
G    = Link.G;
E    = Link.E;
Eta  = Link.Eta;
for ii=1:np
    Ms((ii-1)*6+1:ii*6,:) = Rho*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
    Es((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
    Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
end
   


%eof