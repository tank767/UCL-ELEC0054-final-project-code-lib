%Constants and variables
h = 6.6261*10^-34;
c = 3*10^8;
lambda = 1.55*10^-6;
f = c/lambda;
nsp = 1.58;
G = 16;
ns = 240;
gamma = 1.2*10^-3;
L = 80;
alpha = 0.0461;
Ptx_dbm = -3;
beta2 = -2.17*10^-26;
Nch = 100;
Rs = 5*10^10;
Brx = 5*10^10;

%Linearizing and calculating Leff
Leff = (1-exp(-alpha*L))/alpha;
linear_G = 10^(G/10);
Ptx = 10^-3*10^(Ptx_dbm/10);

%Calculating PSD of noise
Sase = h*f*nsp*(linear_G-1);
Snli = ((2/3)^3*gamma^2*(Leff*10^3)*Ptx^3)*(log(pi^2*abs(beta2)*(Leff*10^3)*(Nch*Rs)^2)/(pi*abs(beta2)*Rs^3));

%Calculating SNR
SNR = Ptx/((Sase+Snli)*ns*Brx);

% QAM BER Simulation and Theoretical Calculation in AWGN
N = 1e7;                  % Number of bits
M_range = [4 8 16 32];                    % MPSK
SNR_dB = 10*log10(SNR);              % Signal-to-Noise Ratio in dB

SER_simulated = zeros(size(M_range));
SER_theoretical = zeros(size(M_range));
BER_theoretical = zeros(size(M_range));

for M_idx = 1:length(M_range)

    % Generate random bits
    bits = randi([0 M_range(M_idx)-1], N, 1);
    
    % Modulate bits using QPSK
    modulated_signal = pskmod(bits, M_range(M_idx), pi/M_range(M_idx));

    % Pass the modulated signal through an AWGN channel
    received_signal = awgn(modulated_signal, SNR_dB);
    
    % Demodulate the received signal
    received_bits = pskdemod(received_signal, M_range(M_idx), pi/M_range(M_idx));
    
    % Calculate the simulated BER
    SER_simulated(M_idx) = sum(bits ~= received_bits)/N;
    
    % Calculate the theoretical BER for QPSK in AWGN
    BER_theoretical(M_idx) = (1/log2(M_range(M_idx)))*erfc(sqrt(SNR*Brx/Rs)*sin(pi/M_range(M_idx)));
    SER_theoretical(M_idx) = 1-(1-BER_theoretical(M_idx))^log2(M_range(M_idx));

    fprintf('SNR = %d dB, Simulated BER = %e, Theoretical BER = %e\n', ...
            SNR_dB, SER_simulated(M_idx), SER_theoretical(M_idx));
end
fprintf('\nSimulation complete. Plotting BER curve.\n');
figure;
semilogy(M_range, SER_simulated, 'bo-', 'DisplayName', 'Simulated PSK');
hold on;
semilogy(M_range, SER_theoretical, 'rx-', 'DisplayName', 'Theoretical PSK');
xlabel('M number');
ylabel('Symbol Error Rate (SER)');
title('PSK Modulation over AWGN Channel');
grid on;
ylim([1e-4 1]);
legend('show');
