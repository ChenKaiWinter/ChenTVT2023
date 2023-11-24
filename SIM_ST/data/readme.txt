% downlink isac
a_tx_1 = zeros(N_t, length(angleSpace));
a_tx_2 = zeros(N_t, length(angleSpace));
a_tx_tgt = zeros(N_t, length(angleSpace));
for i = 1:N_t
    a_tx_1 = a_tx_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(1)/180)));
    a_tx_2 = a_tx_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(2)/180)));
    a_tx_tgt = a_tx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_tx=(a_tx_1+a_tx_2+a_tx_tgt)/3;

% uplink isac
a_rx_1 = zeros(N_r, length(angleSpace));
a_rx_2 = zeros(N_r, length(angleSpace));
a_rx_3 = zeros(N_r, length(angleSpace));
a_rx_4 = zeros(N_r, length(angleSpace));
a_rx_tgt = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx_1 = a_rx_1 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(1)/180)));
    a_rx_2 = a_rx_2 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(2)/180)));
    a_rx_3 = a_rx_3 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(3)/180)));
    a_rx_4 = a_rx_4 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(4)/180)));
    a_rx_tgt = a_rx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_rx=(a_rx_1+a_rx_2+a_rx_3+a_rx_4+a_rx_tgt)/5;

% uplink comm
a_rx_1 = zeros(N_r, length(angleSpace));
a_rx_2 = zeros(N_r, length(angleSpace));
a_rx_3 = zeros(N_r, length(angleSpace));
a_rx_4 = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx_1 = a_rx_1 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(1)/180)));
    a_rx_2 = a_rx_2 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(2)/180)));
    a_rx_3 = a_rx_3 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(3)/180)));
    a_rx_4 = a_rx_4 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(4)/180)));
end
a_rx=(a_rx_1+a_rx_2+a_rx_3+a_rx_4)/4;

% uplink sensing
a_rx_tgt = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx_tgt = a_rx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_rx=a_rx_tgt;

% conv 
SINR_24=[1.79 0.2e+8 3.96e+8 6.07e+8 6.07e+8 6.07e+8 6.07e+8 6.07e+8 6.07e+8 6.07e+8];
SINR_24dBm=plot(10*log10(SINR_24), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
SINR_20=[0.85 0.2e+7 0.5e+8 2.57e+8 2.57e+8 2.57e+8 2.57e+8 2.57e+8 2.57e+8 2.57e+8];
SINR_20dBm=plot(10*log10(SINR_20), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
SINR_16=[0.62 2.8e+5 4.6e+6 2.8e+7 2.8e+7 2.8e+7 2.8e+7 2.8e+7 2.8e+7 2.8e+7];
SINR_16dBm=plot(10*log10(SINR_16), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
legend([SINR_24dBm, SINR_20dBm, SINR_16dBm],'P_t=24dBm', 'P_t=20dBm', 'P_t=16dBm'); 
xlabel('Iteration');
ylabel('Sensing SINR (dB)');

% tradeoff 
SINR_40=[7.01e+9 4.95e+9 2.89e+9 6.26e+8 3.47e+8];
SINR_40dBm=plot(10*log10(SINR_40), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
SINR_30=[6.15e+9 4.11e+9 2.15e+9 4.63e+8 2.19e+8];
SINR_30dBm=plot(10*log10(SINR_30), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
SINR_20=[5.13e+9 2.82e+9 7.90e+8 2.92e+8 1.49e+8];
SINR_20dBm=plot(10*log10(SINR_20), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
legend([SINR_40dBm, SINR_30dBm, SINR_20dBm],'\eta_{RSI}=-40dBm', '\eta_{RSI}=-30dBm', '\eta_{RSI}=-20dBm'); 
xlabel('\gamma_d');
ylabel('Sensing SINR (dB)');

% nomaoma: 
uue1noma=[14.80 13.00 12.40 10.80 10.00];
uue1nomadBm=plot(10*log10(uue1noma), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o', 'Color', 'b');
hold on;
uue2noma=[12.80 11.80 10.20 9.60 8.00];
uue2nomadBm=plot(10*log10(uue2noma),  'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 's', 'Color', 'b');
hold on;
uue3noma=[11.00 9.80 8.20 6.80 6.00];
uue3nomadBm=plot(10*log10(uue3noma),  'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'd', 'Color', 'b');
hold on;
uue4noma=[8.90 7.20 6.70 5.10 4.00];
uue4nomadBm=plot(10*log10(uue4noma),  'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^', 'Color', 'b');
hold on;
uue1oma=[14.80 13.00 12.40 10.80 10.00];
uue1omadBm=plot(10*log10(uue1oma), 'LineWidth', 1.0, 'Linestyle', '--', 'Marker', 'o', 'Color', 'r');
hold on;
uue2oma=[1.96 1.00  0.48 0.30 0.12];
uue2omadBm=plot(10*log10(uue2oma), 'LineWidth', 1.0, 'Linestyle', '--', 'Marker', 's', 'Color', 'r');
hold on;
uue3oma=[0.74 0.4 0.22 0.12 0.04];
uue3omadBm=plot(10*log10(uue3oma), 'LineWidth', 1.0, 'Linestyle', '--', 'Marker', 'd', 'Color', 'r');
hold on;
uue4oma=[0.38 0.18 0.12 0.06 0.012];
uue4omadBm=plot(10*log10(uue4oma), 'LineWidth', 1.0, 'Linestyle', '--', 'Marker', '^', 'Color', 'r');
hold on;
grid on;
xlabel('Transmit Power P_t (dBm)');
ylabel('Uplink Communication SINR (dB)');

% FD
% downlink sinr
50.81 2983.53
% uplink sinr
0.93 0.13 0.0017 0.00037                      
% HD
% downlink sinr
784.20 924.33
% uplink sinr
62.63 58.57  23.56 27.24
% plot  
% RSI\sum rate
% -10 326.64 95.29   10.80 9.60   6.80   5.10   
% -20 449.22 62.22   12.40 10.20 8.20   6.70   
% -30 532.01 362.02 13.91 15.90 13.0   7.60   
% -40 554.82 488.90 23.35 16.66 14.76 8.15   
% -50 599.26 767.09 40.35 18.76 16.32 10.28 
sumRate10=  [326.64 95.29   10.80 9.60   6.80   5.10]; % 27.48
sumRate20=  [449.22 62.22   12.40 10.20 8.20   6.70]; % 28.17
sumRate30=  [532.01 362.02 13.91 15.90 13.0   7.60]; % 32.45
sumRate40=  [554.82 488.90 23.35 16.66 14.76 8.15]; % 33.98
sumRate50=[599.26 767.09 40.35 18.76 16.32 10.28]; % 36.10
sumRate=   [784.20 924.33 62.63 58.57  27.24 23.56]; % 40.80
sum=0;
for i=1:6
sum=sum+log2(1+sumRate(i));
end
FD=[36.10 33.98 32.45 28.17 27.48];
HD=[20.40 20.40 20.40 20.40 20.40];
FDPlot=plot(FD, 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
HDPlot=plot(HD, 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
grid on;
xlabel('RSI (dB)');
ylabel('Sum Rate (bps/Hz)');

% FD_HD_Sensing
HD_SINR=[1.81e+04 1.81e+04 1.81e+04 1.81e+04 1.81e+04];
FD_SINR=[2.28e+9 3.47e+8 2.19e+8 1.49e+8 3.89e+7]
plot(10*log10(HD_SINR), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
plot(10*log10(FD_SINR), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
grid on;
xlabel('RSI (dB)');
ylabel('Sensing SINR (dB)');