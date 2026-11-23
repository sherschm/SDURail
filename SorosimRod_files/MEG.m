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

    J  = zeros(np,1);
    A_p  = zeros(np,1);
    Iy_p = zeros(np,1);
    Iz_p = zeros(np,1);
    Ix_p = zeros(np,1);

    esc_p = zeros(np,1); % shear offset

    for ii = 1:np

        h  = h_fn(Xs(ii));
        w  = w_fn(Xs(ii));

        hw = h - 2*tf;
        l_lip = tw;
        % Areas
        A_web = tw*hw;
        A_fl  = w*tf;
        A_lip = 2*tf*l_lip;

        A = A_web + 2*A_fl + A_lip;

        % Centroid (z-direction)
        z_ref = w/2;
        y_ref = 0;

        % z_web = tw/2;
        % z_fl  = tw + (w-tw)/2;
        % z_lip = w - tf/2;      % lip centroid in z-direction
        z_web = tw/2        - z_ref;
        z_fl  = tw + (w-tw)/2 - z_ref;
        z_lip = w - tf/2    - z_ref;

        zc = ( ...
            A_web*z_web ...
            + 2*A_fl*z_fl ...
            + A_lip*z_lip ) / A;

        % Centroidal Iy
        Iy_web_local = tw*hw^3/12;

        Iy_fl_local = w*tf^3/12;

        y_web = 0;
        y_fl = h/2 - tf/2;
        y_lip = h/2 - tf - l_lip/2;

        Iy_lip_local = tf*l_lip^3/12;

        Iy = ...
            Iy_web_local ...
            + 2*(Iy_fl_local + A_fl*y_fl^2) ...
            + 2*(Iy_lip_local + (tf*l_lip)*y_lip^2);
        % Centroidal Iz
        Iz_web_local = hw*tw^3/12;
        Iz_fl_local  = tf*w^3/12;
        Iz_lip_local = l_lip*tf^3/12;

        % Iz = ...
        %     (Iz_web_local + A_web*(z_web-zc)^2) ...
        %     + 2*(Iz_fl_local + A_fl*(z_fl-zc)^2) ...
        %     + 2*(Iz_lip_local + (tf*l_lip)*(z_lip - zc)^2);

        Iz = ...
            (Iz_web_local + A_web*z_web^2) ...
            + 2*(Iz_fl_local + A_fl*z_fl^2) ...
            + 2*(Iz_lip_local + (tf*l_lip)*z_lip^2);

        % Saint-Venant torsion constant
        Jsv = (1/3)*( ...
            2*w*tf^3 ...
            + hw*tw^3 );

        % Torsional constant computed via Roarks Stress and Strain table 10.3 (1g) 
        t = Link.tw; 
        h_roark = Link.h{1}(0)-t;
        b_roark = Link.w{1}(0)-t/2;

        % % FOR CHANNEL (U) Compute K (Table 10.1, 25 (and 24))
        % a_K = (Link.h{1}(0)/2);
        % c_K = (Link.w{1}(0))-t;
        % K_L1 = a_K*t^3*((1/3)-0.21*(t/a_K)*(1-t^4/(12*(a_K^4))));
        % K_L2 = c_K*t^3*((1/3)-0.105*(t/c_K)*(1-t^4/(192*(c_K^4))));
        % % alpha_K = 0.07;
        % % D_K = 0;
        % K_L = K_L1+K_L2;
        % K = 2*K_L;
        % G = Link.G; %CONFIRM!
        % C_w = (h_roark^2*b_roark^3*t/12)*((2*h_roark+3*b_roark)/(h_roark+6*b_roark));
        % E =Link.E;
        % Beta = sqrt(K*G/(C_w*E));
        % l = Link.L;

        % FOR C-shape Compute K (Table 10.1, 25 (and 24))
        b1_roark = l_lip+t/2;
        e_roark = b_roark*(3*h_roark^2*b_roark + 6*h_roark^2*b1_roark-8*b1_roark^3)/(h_roark^3+6*h_roark^2*b_roark + 6*h_roark^2*b1_roark + 8*b1_roark^3 - 12*h_roark* b1_roark^2);
        K = (t^3/3)*(h_roark+2*b_roark+2*b1_roark);
        G = Link.G; %CONFIRM!
        C_w = t*((h_roark^2*b_roark^2/2)*(b1_roark + b_roark/3 - e_roark - 2*e_roark*b1_roark/b_roark + 2*b1_roark^2/h_roark) ...
            + (h_roark^2*e_roark^2/2)*(b_roark+b1_roark+h_roark/6 - 2*b1_roark^2/h_roark)+(2*b1_roark^3/3)*(b_roark+e_roark)^2);
        E =Link.E;
        Beta = sqrt(K*G/(C_w*E));
        l = Link.L;

        %J(ii) = l*E*C_w*Beta^3/(4*G*((Beta*l)/4-tanh((Beta*l)/4)));

        %Note - Warping absolutely dominates, since Beta*l is very small!
        J(ii) = 48*C_w/l^2;

        A_p(ii)  = A;
        Iy_p(ii) = Iy;
        Iz_p(ii) = Iz;
        Ix_p(ii) = J(ii);

        %For U Channel
        % esc = ...
        %     3*w*tf*(h-tf)^2 / ...
        %     (6*w*tf*(h-tf) + hw*tw^2);

        esc = ...
            ( 3*w*tf*(h-tf)^2 ...
            + 6*l_lip*tf*(w-tf/2)*(h-tf) ) ...
            /( 6*w*tf*(h-tf) ...
            + 2*l_lip*tf*(h-tf) ...
            + hw*tw^2 );

        %esc_p(ii) = esc;
        %esc_p(ii) = 0.0;

    end
end

Ms = zeros(6*np,6); %inertia
Es = zeros(6*np,6); %stiffness
Gs = zeros(6*np,6); %damping

Rho  = Link.Rho;
G    = Link.G;
E    = Link.E;
Eta  = Link.Eta;
if Link.CS =='SDU'   % Channel ("C") section

    Jm_p = Iy_p + Iz_p;   % mass polar inertia
    Js_p = J;             % Saint-Venant torsion constant
    for ii=1:np
        Ms((ii-1)*6+1:ii*6,:) = Rho*diag([Jm_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
        %Es((ii-1)*6+1:ii*6,:) = diag([G*Js_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
       % e = esc_p(ii);

        K = zeros(6);

        K(1,1) = G*Js_p(ii);%+ G*A_p(ii)*e^2;

        K(2,2) = E*Iy_p(ii);
        K(3,3) = E*Iz_p(ii);

        K(4,4) = E*A_p(ii);

        %kappa = 0.3; %Timoshenko shear
        kappa = 0.5; %Timoshenko shear

        %kappa = 1.0; %Timoshenko shear

        K(5,5) = kappa*G*A_p(ii);
        K(6,6) = kappa*G*A_p(ii);

        %K(1,6) = -G*A_p(ii)*e;
       % K(6,1) = -G*A_p(ii)*e;

        Es((ii-1)*6+1:ii*6,:) = K;
        Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Jm_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
    end

else
    for ii=1:np
        Ms((ii-1)*6+1:ii*6,:) = Rho*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
        Es((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
        Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
    end
end
% 
% %This function computes Ms, Es, and Gs for a given Xs for the jth division
% %of a link (15.12.2022)
% 
% function [Ms,Es,Gs]= MEG(Link,j,Xs)
% np = length(Xs);
% %j is division number
% 
% if Link.CS=='R'
%     h_fn     = Link.h{j};
%     w_fn     = Link.w{j};
% 
%     %updating:
%     h_nGauss = zeros(np,1);
%     w_nGauss = zeros(np,1);
%     for ii=1:np
%         h_nGauss(ii) = h_fn(Xs(ii));
%         w_nGauss(ii) = w_fn(Xs(ii));
%     end
% 
%     Iy_p = (1/12)*h_nGauss.*(w_nGauss.^3);
%     Iz_p = (1/12)*(h_nGauss.^3).*w_nGauss;
%     Ix_p = Iy_p+Iz_p;
%     A_p  = h_nGauss.*w_nGauss;
% 
% elseif Link.CS=='C' %Circular
%     r_fn     = Link.r{j};
% 
%     %updating:
%     r_nGauss = zeros(np,1);
%     for ii=1:np
%         r_nGauss(ii) = r_fn(Xs(ii));
%     end
%     Iy_p = (pi/4)*r_nGauss.^4;
%     Iz_p = Iy_p;
%     Ix_p = Iy_p+Iz_p;
%     A_p  = pi*r_nGauss.^2;
% 
% elseif Link.CS=='E'
%     a_fn     = Link.a{j};
%     b_fn     = Link.b{j};
% 
%     %updating:
%     a_nGauss = zeros(np,1);
%     b_nGauss = zeros(np,1);
%     for ii=1:np
%         a_nGauss(ii) = a_fn(Xs(ii));
%         b_nGauss(ii) = b_fn(Xs(ii));
%     end
% 
%     Iy_p = (pi/4)*a_nGauss.*(b_nGauss.^3);
%     Iz_p = (pi/4)*(a_nGauss.^3).*b_nGauss;
%     Ix_p = Iy_p+Iz_p;
%     A_p  = pi*a_nGauss.*b_nGauss;
% 
% elseif Link.CS=='SDU'   % Channel ("C") section
% 
%     h_fn  = Link.h{j};
%     w_fn  = Link.w{j};
%     tw = Link.tw;
%     tf = Link.tf;
% 
%     A_p  = zeros(np,1);
%     Iy_p = zeros(np,1);
%     Iz_p = zeros(np,1);
%     Ix_p = zeros(np,1);
% 
%     for ii = 1:np
% 
%         h  = h_fn(Xs(ii));
%         w  = w_fn(Xs(ii));
% 
%         hw = h - 2*tf;
% 
%         % Areas
%         A_web = tw*hw;
%         A_fl  = w*tf;
% 
%         A = A_web + 2*A_fl;
% 
%         % Centroid (z-direction)
%         z_web = tw/2;
%         z_fl  = tw + (w-tw)/2;
% 
%         zc = ( ...
%             A_web*z_web ...
%             + 2*A_fl*z_fl ) / A;
% 
%         % Centroidal Iy
%         Iy_web = tw*hw^3/12;
% 
%         Iy_fl_local = w*tf^3/12;
%         y_fl = h/2 - tf/2;
% 
%         Iy = Iy_web ...
%             + 2*(Iy_fl_local + A_fl*y_fl^2);
% 
%         % Centroidal Iz
%         Iz_web_local = hw*tw^3/12;
%         Iz_fl_local  = tf*w^3/12;
% 
%         Iz = ...
%             (Iz_web_local + A_web*(z_web-zc)^2) ...
%             + 2*(Iz_fl_local + A_fl*(z_fl-zc)^2);
% 
%         % Saint-Venant torsion constant
%         J = (1/3)*( ...
%             2*w*tf^3 ...
%             + hw*tw^3 );
% 
%         A_p(ii)  = A;
%         Iy_p(ii) = Iy;
%         Iz_p(ii) = Iz;
%         Ix_p(ii) = J;
%     end
% 
% end
% 
% Ms = zeros(6*np,6); %inertia
% Es = zeros(6*np,6); %stiffness
% Gs = zeros(6*np,6); %damping
% 
% Rho  = Link.Rho;
% G    = Link.G;
% E    = Link.E;
% Eta  = Link.Eta;
% for ii=1:np
%     Ms((ii-1)*6+1:ii*6,:) = Rho*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
%     Es((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
%     Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
% end
% 
% 
% 
% %eof