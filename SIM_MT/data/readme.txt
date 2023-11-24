% downlink isac
a_tx_1 = zeros(N_t, length(angleSpace));
a_tx_2 = zeros(N_t, length(angleSpace));
a_tx_tgt_1 = zeros(N_t, length(angleSpace));
a_tx_tgt_2 = zeros(N_t, length(angleSpace));
a_tx_tgt_3 = zeros(N_t, length(angleSpace));
for i = 1:N_t
    a_tx_1 = a_tx_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(1)/180)));
    a_tx_2 = a_tx_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(2)/180)));
    a_tx_tgt_1 = a_tx_tgt_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(1)/180)));
    a_tx_tgt_2 = a_tx_tgt_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(2)/180)));
    a_tx_tgt_3 = a_tx_tgt_3 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(3)/180)));
end
a_tx=(a_tx_1+a_tx_2+a_tx_tgt_1+a_tx_tgt_2+a_tx_tgt_3)/5;

% uplink isac
a_rx_1 = zeros(N_r, length(angleSpace));
a_rx_2 = zeros(N_r, length(angleSpace));
a_rx_3 = zeros(N_r, length(angleSpace));
a_rx_4 = zeros(N_r, length(angleSpace));
a_rx_tgt_1 = zeros(N_r, length(angleSpace));
a_rx_tgt_2 = zeros(N_r, length(angleSpace));
a_rx_tgt_3 = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx_1 = a_rx_1 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(1)/180)));
    a_rx_2 = a_rx_2 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(2)/180)));
    a_rx_3 = a_rx_3 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(3)/180)));
    a_rx_4 = a_rx_4 + exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*UUE_agl(4)/180)));
    a_rx_tgt_1 = a_rx_tgt_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(1)/180)));
    a_rx_tgt_2 = a_rx_tgt_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(2)/180)));
    a_rx_tgt_3 = a_rx_tgt_3 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(3)/180)));
end
a_rx=(a_rx_1+a_rx_2+a_rx_3+a_rx_4+a_rx_tgt_1+a_rx_tgt_2+a_rx_tgt_3)/7;

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
a_rx_tgt_1 = zeros(N_r, length(angleSpace));
a_rx_tgt_2 = zeros(N_r, length(angleSpace));
a_rx_tgt_3 = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx_tgt_1 = a_rx_tgt_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(1)/180)));
    a_rx_tgt_2 = a_rx_tgt_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(2)/180)));
    a_rx_tgt_3 = a_rx_tgt_3 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(3)/180)));
end
a_rx=(a_rx_tgt_1+a_rx_tgt_2+a_tgt_3)/3;

% conv
figure
zeta_24=[7500 6250 6875.00 7187.50 7031.250 6953.1250 6992.18750 7011.718750 7001.9531250 6997.07031250 6994.628906250 6993.4082031250 6992.79785156250 6992.49267578125 6992.64526367188 6992.72155761719 6992.72155761719 6992.72155761719 6992.72155761719 6992.72155761719];
zeta_20=[7500 6250 5625.00 5312.50 5156.250 5078.1250 5039.06250 5019.531250 5009.7656250 5004.88281250 5002.441406250 5001.2207031250 5000.61035156250 5000.30517578125 5000.15258789063 5000.07629394531 5000.03814697266 5000.03814697266 5000.03814697266 5000.03814697266];
zeta_16=[2500 1875 2187.50 2343.75 2265.625 2304.6875 2324.21875 2333.984375 2329.1015625 2331.54296875 2330.322265625 2330.9326171875 2331.23779296875 2331.39038085938 2331.31408691406 2331.27593994141 2331.25686645508 2331.24732971191 2331.24256134033 2331.24017715454];
zeta_24dBm=plot(10*log10(zeta_24),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
zeta_20dBm=plot(10*log10(zeta_20),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 's');
hold on;
grid on;
zeta_16dBm=plot(10*log10(zeta_16),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'd');
hold on;
grid on;
xlabel('Iteration');
ylabel('Sensing SINR (dB)');
legend([zeta_24dBm, zeta_20dBm, zeta_16dBm],'P_t=24dBm', 'P_t=20dBm', 'P_t=16dBm');

% tradeoff
figure;
SINR_40=[6.46e+06 4.39e+06 3.34e+06 3.13e+06 2.35e+06];
SINR_30=[1.99e+06 1.73e+06 1.41e+06 1.30e+06 1.24e+06];
SINR_20=[1.21e+06 1.00e+06 9.08e+05 7.03e+05 6.09e+05];
SINR_40dBm=plot(10*log10(SINR_40),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
SINR_30dBm=plot(10*log10(SINR_30),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 's');
hold on;
grid on;
SINR_20dBm=plot(10*log10(SINR_20),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'd');
hold on;
grid on;
xlabel('\gamma_{k,d}');
ylabel('Sensing SINR (dB)');
legend([SINR_40dBm, SINR_30dBm, SINR_20dBm],'\sigma _{RSI}^2=-40dBm', '\sigma _{RSI}^2=-30dBm', '\sigma _{RSI}^2=-20dBm');

% nomaoma
% noma
P_H=zeros(N_r,N_r,M);
for tilde_m=1:M
    if tilde_m==M
        P_H(:,:,tilde_m)=zeros(N_r,N_r);
    end
    for m=tilde_m+1:M
        P_H(:,:,tilde_m)=P_H(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
    end
end
SINR_u=zeros(1,4);
for m=1:M
    SINR_u(m)=P_u*abs(real(v_m(:,m)'*h_u(:,m)))^2/real( v_m(:,m)'*(P_H(:,:,m)+R_e+R_n+n_u*eye(N_r))*v_m(:,m) );
end
% oma
P_H_oma=zeros(N_r,N_r,M);
for tilde_m=1:M
    for m=1:M
        if(m~=tilde_m)
            P_H_oma(:,:,tilde_m)=P_H(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
        end
    end
end
SINR_u_oma=zeros(1,4);
for m=1:M
    SINR_u_oma(m)=P_u*abs(real(v_m(:,m)'*h_u(:,m)))^2/real( v_m(:,m)'*(P_H_oma(:,:,m)+R_e+R_n+n_u*eye(N_r))*v_m(:,m) );
end
figure;
SINR_u1_noma=[13.98 12.60 11.40 10.60 10.00];
SINR_u2_noma=[12.40 11.30 10.10 9.20 8.00];
SINR_u3_noma=[10.96 9.50 8.00 6.60 6.00];
SINR_u4_noma=[8.80 6.80 5.70 4.90 4.00];
SINR_u1_oma=[13.98 12.60 11.40 10.60 10.00];
SINR_u2_oma=[1.66e-03 5.40e-04 5.14e-04 3.00e-04 1.57e-04];
SINR_u3_oma=[1.53e-04 1.05e-04 2.73e-05 2.31e-05 1.75e-05];
SINR_u4_oma=[8.01e-06 5.38e-06 1.69e-06 1.51e-06 1.29e-06];
SINR_u1_noma_dB=plot(10*log10(SINR_u1_noma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
SINR_u2_noma_dB=plot(10*log10(SINR_u2_noma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 's');
hold on;
grid on;
SINR_u3_noma_dB=plot(10*log10(SINR_u3_noma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'd');
hold on;
grid on;
SINR_u4_noma_dB=plot(10*log10(SINR_u4_noma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
grid on;
SINR_u1_oma_dB=plot(10*log10(SINR_u1_oma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'o');
hold on;
grid on;
SINR_u2_oma_dB=plot(10*log10(SINR_u2_oma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 's');
hold on;
grid on;
SINR_u3_oma_dB=plot(10*log10(SINR_u3_oma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', 'd');
hold on;
grid on;
SINR_u4_oma_dB=plot(10*log10(SINR_u4_oma),'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');
hold on;
grid on;
xlabel('Transmit Power P_t (dBm)');
ylabel('Uplink Communication SINR (dB)');

% Sensing SINR I&P
FD=[4.79e+03 3.04e+03 1.57e+03 1.17e+03  8.02e+02];
plot(10*log10(FD), 'LineWidth', 1.0, 'Linestyle', '-', 'Marker', '^');