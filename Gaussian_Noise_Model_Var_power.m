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
Ptx_dbm_range = -10:0.5:10;
beta2 = -2.17*10^-26;
Nch = 100;
Rs = 5*10^10;
Brx = 5*10^10;

%Linearizing and calculating Leff
Leff = (1-exp(-alpha*L))/alpha;
linear_G = 10^(G/10);

SER_simulated = zeros(size(Ptx_dbm_range));
SER_theoretical = zeros(size(Ptx_dbm_range));

for ptx_idx = 1:length(Ptx_dbm_range)
    current_ptx = 10^-3*10^(Ptx_dbm_range(ptx_idx)/10);
    %Calculating PSD of noise
    Sase = h*f*nsp*(linear_G-1);
    Snli = ((2/3)^3*gamma^2*(Leff*10^3)*current_ptx^3)*(log(pi^2*abs(beta2)*(Leff*10^3)*(Nch*Rs)^2)/(pi*abs(beta2)*Rs^3));
    %Calculating SNR
    SNR = current_ptx/((Sase+Snli)*ns*Brx);
    
    % QPSK BER Simulation and Theoretical Calculation in AWGN
    
    N = 1e7;                  % Number of bits
    M = 4;
    SNR_dB = 10*log10(SNR);              % Signal-to-Noise Ratio in dB
    
    % Generate random bits
    bits = randi([0 M-1], N, 1);
    
    % Modulate bits using QPSK
    modulated_signal = pskmod(bits, M, pi/M);
    
    % Pass the modulated signal through an AWGN channel
    received_signal = awgn(modulated_signal, SNR_dB);
    
    % Demodulate the received signal
    received_bits = pskdemod(received_signal, M ,pi/M);
    
    % Calculate the simulated BER
    SER_simulated(ptx_idx) = sum(bits ~= received_bits)/N;
    
    % Calculate the theoretical BER for QPSK in AWGN
    BER_theoretical = (1/log2(M))*erfc(sqrt(SNR*Brx/Rs)*sin(pi/M));
    SER_theoretical(ptx_idx) = 1-(1-BER_theoretical)^log2(M);

    fprintf('SNR = %d dB, Simulated BER = %e, Theoretical BER = %e\n', ...
            SNR_dB, SER_simulated(ptx_idx), SER_theoretical(ptx_idx));
end
fprintf('\nSimulation complete. Plotting BER curve.\n');
figure;
semilogy(Ptx_dbm_range, SER_simulated, 'bo-', 'DisplayName', 'Simulated QPSK');
hold on;
semilogy(Ptx_dbm_range, SER_theoretical, 'r--', 'DisplayName', 'Theoretical QPSK');
xlabel('transmisstion power');
ylabel('Symbol Error Rate (SER)');
title('QPSK Modulation over AWGN Channel');
grid on;
ylim([1e-4 1]);
legend('show');
% Display results