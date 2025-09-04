%Constants and variables
h = 6.6261*10^-34;
c = 3*10^8;
lambda = 1.52*10^-6;
f = c/lambda;
nsp = 1.58;
G = 16;
ns = 8;
gamma = 1.2*10^-3;
L = 80;
alpha = 0.0461;
Ptx_dbm = -1.5;
beta2 = -2.85*10^-26;
Nch = 800;
Rs = 5*10^10;
Brx = 5*10^10;

%Linearizing and calculating Leff
Leff = (1-exp(-alpha*L))/alpha;
linear_G = 10^(G/10);
Ptx = 10^-3*10^(Ptx_dbm/10);
Ptx_dbw = Ptx_dbm-30;

%Calculating PSD of noise
Sase = h*f*nsp*(linear_G-1);
Snli = ((2/3)^3*gamma^2*(Leff*10^3)*Ptx^3)*(log(pi^2*abs(beta2)*(Leff*10^3)*(Nch*Rs)^2)/(pi*abs(beta2)*Rs^3));
%Calculating SNR
SNR = Ptx/((Sase+Snli)*ns*Brx);

% QPSK BER Simulation and Theoretical Calculation in AWGN

N = 1e7;                  % Number of bits
M = 64;                    % M-QAM
SNR_dB = 10*log10(SNR);              % Signal-to-Noise Ratio in dB

% Generate random bits
bits = randi([0 M-1], N, 1);

% Modulate bits using QPSK
modulated_signal = qammod(bits, M);

%average power normalization
%minD = helperAvgPow2MinD(Ptx,M);
%tx_signal = (minD/2) .* modulated_signal;

% Pass the modulated signal through an AWGN channel
received_signal = awgn(modulated_signal, SNR_dB,"measured");

% Demodulate the received signal
received_bits = qamdemod(received_signal, M);

% Calculate the simulated BER
SER_simulated = sum(bits ~= received_bits)/N;

% Calculate the theoretical BER for QPSK in AWGN
BER_theoretical = (2/log2(M))*(1-1/sqrt(M))*erfc(sqrt(1.5*SNR/(M-1)));
SER_theoretical = 1-(1-BER_theoretical)^log2(M);
% biterate
Biterate = log2(M)*Rs*Nch;

% Display results
fprintf('Simulated SER: %e\n', SER_simulated);
fprintf('Theoretical SER: %e\n', SER_theoretical);
fprintf('Bitreate: %e\n', Biterate);
