%% Inputs
% carriage_motion : Nt x Ndof matrix
% dt       : sample time [s]

[Nt, Ndof] = size(carriage_motion);

Fs = 1/dt;              % Sampling frequency [Hz]
f = Fs*(0:floor(Nt/2))/Nt;

% FFTs
amplitude = zeros(length(f), Ndof);
dominant_freq = zeros(Ndof,1);
dominant_amp  = zeros(Ndof,1);

idx_free = 0.02/dt;
for j = 1:Ndof

    x = carriage_motion(idx_free:end,j);
    
    % Remove mean (important)
    x = x - mean(x);

    % FFT
    Y = fft(x);

    % Two-sided spectrum
    P2 = abs(Y/Nt);

    % Single-sided spectrum
    P1 = P2(1:floor(Nt/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);

    amplitude(:,j) = P1;

    % Ignore zero-frequency peak
    [dominant_amp(j), idx] = max(P1(2:end));

    idx = idx + 1;
    dominant_freq(j) = f(idx);

    %torsional_freq(j) = 
end

%Bending
I = RailLink.w{1}(0)*RailLink.h{1}(0)^3/12;
EI = RailLink.E*I;
m = RailCarriage.VLinks(1).M(4,4);
f1_bending_analytical = (1/(2*pi))*sqrt(192*EI/(m*RailLink.L^3));

%%
% Torsion
I_roll = RailCarriage.VLinks(1).M(1,1);
%J = RailLink.w{1}(0)*RailLink.h{1}(0)^3/3;

J = 0.1406* RailLink.w{1}(0)^4;
%J = (1/6)* RailLink.w{1}(0)^4;
G = RailLink.G;
f1_torsion_analytical = (1/(2*pi))*sqrt(4*G*J/(I_roll*RailLink.L));

x = carriage_motion(idx_free:end,j);

% Remove mean (important)
x = x - mean(x);

% FFT
Y = fft(x);

% Two-sided spectrum
P2 = abs(Y/Nt);

% Single-sided spectrum
P1 = P2(1:floor(Nt/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

amplitude(:,j) = P1;

% Ignore zero-frequency peak
[dominant_amp(j), idx] = max(P1(2:end));

idx = idx + 1;
dominant_freq(j) = f(idx);


%lateral deformation

g_s0 = FwdKinematicsAtS(RailLinkage,q_equil,s0); %SE(3)
z_carriage_equil = g_s0(3,4);

z_carriage_equil_analytical = m*9.81*RailLink.L^3/(192*EI);

%%

save(fullfile('data\dynamic_data\slender_scenario\shape_functions_vs_nat_freq\','data.mat'), ...
    'f','amplitude','z_carriage_equil','dominant_freq','carriage_motion');

close all
figure('Color','w');

ax = axes;
hold(ax,'on')

% White axes background
set(ax,...
    'Color','w',...
    'XColor','k',...
    'YColor','k',...
    'FontName','Times',...
    'FontSize',11,...
    'LineWidth',1.0);

grid on
box on

% Responses
plot(f, amplitude(:,1)./max(abs(amplitude(:,1))), ...
    'LineWidth',1.3,'Color', [0.8500 0.3250 0.0980]);

plot(f, amplitude(:,6)./max(abs(amplitude(:,6))), ...
    'LineWidth',1.3,'Color', [0 0.4470 0.7410]);

% Analytical frequencies
xline(f1_torsion_analytical,'--','LineWidth',1.8,'Color', [0.8500 0.3250 0.0980]);
xline(f1_bending_analytical,'--','LineWidth',1.8,'Color', [0 0.4470 0.7410]);

xlabel('Frequency (Hz)', ...
    'FontSize',12, ...
    'Color','k');

ylabel('Normalised Amplitude (-)', ...
    'FontSize',12, ...
    'Color','k');

legend({'Roll response',...
    'Vertical response',...
    'Analytical torsion',...
    'Analytical bending'}, ...
    'Location','best');

set(legend,'Color','w')
set(legend,'TextColor','k')
set(legend,'EdgeColor','k')

% Squish vertically
ylim([0 1.1])      % choose range of interest
xlim([0 150])      % choose range of interest

% ylim([0 0.25])   % even flatter
% ylim([0 0.15])   % highlights small peaks

% Make figure itself shorter
set(gcf,'Position',[100 100 700 250]);

set(gcf,'Renderer','painters');

exportgraphics(gcf, ...
    'plots/carriage_FFT.pdf', ...
    'ContentType','vector');