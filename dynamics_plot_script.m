%% Inputs
% carriage_motion : Nt x Ndof matrix
% dt       : sample time [s]

[Nt, Ndof] = size(carriage_motion);

Fs = 1/dt;              % Sampling frequency [Hz]
f = Fs*(0:floor(Nt/2))/Nt;

%% FFTs
amplitude = zeros(length(f), Ndof);
dominant_freq = zeros(Ndof,1);
dominant_amp  = zeros(Ndof,1);

for j = 1:Ndof

    x = carriage_motion(:,j);

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

end

%% Plot all DOFs
figure;
hold on

for j = 1:Ndof
    plot(f, amplitude(:,j), 'LineWidth', 1);
end

xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('FFT Spectrum')
grid on

%% Report dominant frequencies
fprintf('\nDominant frequencies:\n');
for j = 1:Ndof
    fprintf('DOF %2d : %8.3f Hz\n', j, dominant_freq(j));
end

ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontSize = 11;
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;


addpath(".\data\dynamic_data\")

T = readmatrix('lateral_response.xlsx');
plot(T(1,2:end),T(30,2:end))
%%
T2 = readmatrix('payload_a.csv');

% Original data
t = T2(2:end,2);
y = T2(2:end,3)-mean(T2(2:end,3));

% Uniform time vector
dt = 0.01;
t_uniform = t(1):dt:t(end);

% Interpolate response onto uniform grid
y_uniform = interp1(t, y, t_uniform, 'pchip');

% Plot comparison
figure
%plot(t, y, 'o', t_uniform, y_uniform, '-')
plot(t_uniform, (y_uniform-y_uniform(1))*1e-3, '-')

hold on 
plot(tvec_out, carriage_motion, '-')
xlabel('Time (s)')
legend('ANSYS','Cosserat')
ylabel('Response')
grid on

Nt = length(y_uniform)-1;

Fs = 1/dt;              % Sampling frequency [Hz]
f = Fs*(0:floor(Nt/2))/Nt;
amplitude = zeros(length(f), Ndof);
for j = 1:Ndof
    % ansys
    x = T2(2:end,3);

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

end