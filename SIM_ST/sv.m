function [a_rx,a_tx] = sv(N_r,N_t,lambda,angleSpace,TxAntLoc,RxAntLoc,DUE_agl,UUE_agl,tgt_agl,tgt_index)

a_tx_1 = zeros(N_t, length(angleSpace));
a_tx_2 = zeros(N_t, length(angleSpace));
a_tx_tgt = zeros(N_t, length(angleSpace));
for i = 1:N_t
    a_tx_1 = a_tx_1 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(1)/180)));
    a_tx_2 = a_tx_2 + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*DUE_agl(2)/180)));
    a_tx_tgt = a_tx_tgt + exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)-sin(pi*tgt_agl(tgt_index)/180)));
end
a_tx=(a_tx_1+a_tx_2+a_tx_tgt)/3;

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
end