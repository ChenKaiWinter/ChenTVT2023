% "Signal Processing and Beamforming Optimization for Full-Duplex NOMA Enabled Integrated Sensing and Communication" 
% Kai Chen, 09/23/2023

%% Set Up the System -- Singel Target
K = 2; DUE_dis = [128 128]; DUE_agl = [-50 50];           % K: number of DUE
M = 4; UUE_dis = [8 16 32 64]; UUE_agl = [-70 -30 40 60]; % M: number of UUE
L = 3; tgt_agl = [-10 0 20]; tgt_index=2;                 % L: number of tgt

N_t = 16;  % tranmit antennas number
N_r = 16;  % receive antennas number

P_t_dB = 20; P_t = 10^(P_t_dB/10); % transmit power at BS
P_u_dB = 10; P_u = 10^(P_u_dB/10); % transmit power at UUE

n_d_dB = -80; n_d = 10^(n_d_dB/10);            % downlink noise power sigma_d
n_u_dB = -80; n_u = 10^(n_u_dB/10);            % uplink noise power sigma_u
n_dB = -80; n = 10^(n_dB/10);                  % uplink noise power sigma after vectorization
n_l = [10^(0/10) 10^(0/10) 10^(0/10)];         % complex amplitude variance |alpha_l|^2
% n_l = [10^(-30/10) 10^(-30/10) 10^(-30/10)]; % complex amplitude variance |alpha_l|^2

f_c = 5e9; c = 3e8; lambda = c/f_c; spacing = lambda/2; % half wavelength
TxAntLoc = spacing*[0:N_t-1];                           % transmit antennas spatial position
RxAntLoc = spacing*[0:N_r-1];                           % receive antennas spatial position
D = 181; angleSpace = linspace(-pi/2, pi/2, D); angleSpaceDeg = linspace(-90, 90, D);

% transmit antennas steering vecter all angle
a_tx = zeros(N_t, length(angleSpace));
for i = 1:N_t
    a_tx(i,:) = exp((1j * 2 * pi * TxAntLoc(i) / lambda) .* (sin(angleSpace)));
end

% receive antennas steering vector all angle
a_rx = zeros(N_r, length(angleSpace));
for i=1:N_r
    a_rx(i,:) = exp((1j * 2 * pi * RxAntLoc(i) / lambda) .* (sin(angleSpace)));
end

% guide steering vector
b = zeros(N_t,L); % transmit steering vecter
for l = 1:L
    b(:,l) = a_tx(:,90+tgt_agl(l));
end
a = zeros(N_r,L); % receive steering vecter
for l = 1:L
    a(:,l) = a_rx(:,90+tgt_agl(l));
end

% generate the user's coordinates
DUE_loc_x = zeros(1, K);
DUE_loc_y = zeros(1, K);
UUE_loc_x = zeros(1, M);
UUE_loc_y = zeros(1, M);
for k =1:K
    DUE_loc_x(k) = DUE_dis(k)*cos(DUE_agl(k));
end
for k =1:K
    DUE_loc_y(k) = DUE_dis(k)*sin(DUE_agl(k));
end
DUE_loc = [DUE_loc_x;DUE_loc_y]; % location of DUE
for m =1:M
    UUE_loc_x(m) = UUE_dis(m)*cos(UUE_agl(m));
end
for m =1:M
    UUE_loc_y(m) = UUE_dis(m)*sin(UUE_agl(m));
end
UUE_loc = [UUE_loc_x;UUE_loc_y]; % location of UUE

% distance between UUE and DUE
D_U_dis = zeros(M,K);
for m = 1:M
    for k = 1:K
        D_U_dis(m,k) = gen_dis(UUE_loc(:,m),DUE_loc(:,k));
    end
end

% w_intial
w_0 = sqrt(P_t/(2*N_t^2))*randn(N_t^2, 1)+sqrt(P_t/(2*N_t^2))*1i*randn(N_t^2, 1);
save('mat\w_0.mat','w_0');
load('mat\w_0.mat');                              % get the inital power

% ch_intial
G_dB = 10;                                        % Rician factor
beta = 3.8;
[a_rx, a_tx]=sv(N_t,N_r,lambda,angleSpace,TxAntLoc,RxAntLoc,DUE_agl,UUE_agl,tgt_agl,tgt_index);
[h_d, h_u, g, H_RSI] = gen_h_rci(G_dB, beta, a_tx, a_rx, N_t, N_r, K, M, ...
    DUE_agl, UUE_agl, DUE_dis, UUE_dis, D_U_dis); % h_d, h_u, g are chanel matrix
save('mat\ch.mat','h_d','h_u','g','H_RSI');       % save the channel
load('mat\ch.mat');                               % get the inital chanel

% generate the matrix that I need
w_k_1=w_0; % w^(k-1)
[C, J_i, R_e, R_n, tilde_R_e, tilde_R_n] = gen_mat(N_r, N_t, L, H_RSI, a, b, n_l, tgt_index, w_k_1);

%% begin
e = 0.001; iter = 1; maxSINR = 0;

while true

% Iter I - u_l: the sensing beamforming vector
% noma
u_l=( (tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1 )...
    /(w_k_1'*kron(eye(N_t),C(:,:,tgt_index))'*(tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1);

gamma_l_k_1=(n_l(tgt_index)*(abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_k_1)^2))/...
    real(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);

% oma
% e & sum{h_u}
% e=ones(N_t,1);
% sum_P_H=zeros(N_r*N_t,N_r*N_t);
% for m=1:M
%     sum_P_H=sum_P_H+kron(h_u(:,m)*h_u(:,m)',e*e');
% end
% u_l=( (P_u*sum_P_H+tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1 )...
%     /(w_k_1'*kron(eye(N_t),C(:,:,tgt_index))'*(P_u*sum_P_H+tilde_R_e+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1);
% 
% gamma_l_k_1=(n_l(tgt_index)*(abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_k_1)^2))/...
%     real(u_l'*(P_u*sum_P_H+tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);

% Iter II - v_m: the uplink receive communication beamforming vector
% noma
P_u_H_u=zeros(N_r,N_r,M);
for tilde_m=1:M
    if tilde_m==M
        P_u_H_u(:,:,tilde_m)=zeros(N_r,N_r);
    end
    for m=tilde_m+1:M
        P_u_H_u(:,:,tilde_m)=P_u_H_u(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
    end
end
v_m=zeros(N_r,M);
for m=1:M
    v_m(:,m)=sqrt(P_u)*(P_u_H_u(:,:,m)+R_e+R_n+n_u*eye(N_r))^(-1)*h_u(:,m);
end

% oma
% P_u_H_u=zeros(N_r,N_r,M);
% for tilde_m=1:M
%     for m=1:M
%         if m~=tilde_m
%             P_u_H_u(:,:,tilde_m)=P_u_H_u(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
%         end
%     end
% end
% v_m=zeros(N_r,M);
% for m=1:M
%     v_m(:,m)=sqrt(P_u)*(P_u_H_u(:,:,m)+R_e+R_n+n_u*eye(N_r))^(-1)*h_u(:,m);
% end

% Iter III - w: the transmit beamforming vector
% O
O=n_l(tgt_index)*(kron(eye(N_t),C(:,:,tgt_index))'*(u_l*u_l')*kron(eye(N_t),C(:,:,tgt_index)));

% Q
n_l_temp=n_l;
n_l_temp(tgt_index)=[];
C_temp=C;
C_temp(:,:,tgt_index)=[];
Q_part_1=zeros(N_r*N_t,N_r*N_t);
for l=1:L-1
    Q_part_1=Q_part_1+n_l_temp(l)*kron(eye(N_t),C_temp(:,:,l))'*(u_l*u_l')*kron(eye(N_t),C_temp(:,:,l));
end
Q=Q_part_1+kron(eye(N_t),H_RSI)'*(u_l*u_l')*kron(eye(N_t),H_RSI)+n*(u_l'*u_l)/P_t;

% tilde_H_k
tilde_H_k_part_2=zeros(N_t^2,N_t^2,K);
for k=1:K
    h_d_temp=h_d;
    h_d_temp(:,k)=[];
    for k_temp=1:K-1
        tilde_H_k_part_2(:,:,k)=tilde_H_k_part_2(:,:,k)+J_i(k_temp)'*h_d_temp(:,k_temp)*h_d_temp(:,k_temp)'*J_i(k_temp);
    end
end
gamma_k_d=[18 18];
tilde_H_k=zeros(N_t^2,N_t^2,K);
for k=1:K
    tilde_H_k(:,:,k)=(J_i(k)'*h_d(:,k)*h_d(:,k)'*J_i(k))/(gamma_k_d(k))-tilde_H_k_part_2(:,:,k);
end

% c_k_d
c_k_d=zeros(1,K);
for k=1:K
    c_k_d(:,k)=P_u*norm(g(:,k))^2+n_d;
end

% P_m
P_m_part_1=zeros(N_t^2,N_t^2,M);
for m=1:M
    for i=1:N_t
        P_m_part_1(:,:,m)=J_i(i)'*H_RSI'*v_m(:,m)*v_m(:,m)'*H_RSI*J_i(i);
    end
end
P_m_part_2=zeros(N_t^2,N_t^2,M);
for m=1:M
    for i=1:N_t
        for l=1:L
            P_m_part_2(:,:,M)=n_l(l).*J_i(i)'*C(:,:,l)'*v_m(:,m)*v_m(:,m)'*C(:,:,l)*J_i(i);
        end
    end
end
P_m=zeros(N_t^2,N_t^2,M);
for m=1:M
    P_m(:,:,m)=P_m_part_1(:,:,m)+P_m_part_2(:,:,m);
end

% c_m_u
c_m_u=zeros(1,m);
gamma_m_u=[10 8 6 4]; 
P_u_V_H=zeros(1,M);
for tilde_m=1:M
    if tilde_m==M
        P_u_V_H(:,tilde_m)=0;
    end
    for m=tilde_m+1:M
        P_u_V_H(:,tilde_m)=P_u_V_H(:,tilde_m)+P_u*abs(v_m(:,m)'*h_u(:,tilde_m))^2;
    end
end
for m=1:M
    c_m_u(:,m)=(P_u*abs(v_m(:,m)'*h_u(:,m))^2)/gamma_m_u(m)-P_u_V_H(:,m)-n_u*(v_m(:,m)'*v_m(:,m));
end

%%%%%%Original Scheme%%%%%% 
cvx_begin SDP quiet
cvx_precision high
variable Z(N_t^2, N_t^2) hermitian
variable y
maximize real(trace(O*Z))
subject to
real(trace(Q*Z))==1;
for k=1:K
    real(trace(tilde_H_k(:,:,k)*Z))>=y*c_k_d(k);
end
for m=1:M
    real(trace(P_m(:,:,m)*Z))<=y*c_m_u(m);
end
trace(Z)<=y*P_t;
Z == semidefinite(N_t^2);
y>=0;
cvx_end

% % eigenvalue decomposition: sqrt(dominant eigenvalue)*dominant eigenvector (It has been proved that there is no rank one solution)
% [V,D] = eigs(Z);
% d = diag(D);
% [mv,mi]= max(d);
% w_dom = sqrt(d(mi))*V(:,mi);
% w_opt = w_dom;

% gaussian randomization -- do randomization to get our feasible beamforming vectors w_cand
tilde_W=Z/y;
nRand = 100;
w_rand = zeros(N_t^2, nRand);

% generate and scale to meet power constraints and other constraints
for L_rand = 1:nRand
    zeta(:,L_rand) = mvnrnd(zeros(N_t^2,1),tilde_W) + 1i*mvnrnd(zeros(N_t^2,1),tilde_W);         % generate nRand
    w_rand(:,L_rand) = sqrt(trace(tilde_W))*zeta(:,L_rand)/sqrt(zeta(:,L_rand)'*zeta(:,L_rand)); % scale them so it adheres to power constraint
end

for i=1:nRand
    if ( ((n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l)>maxSINR) )
        maxSINR=(n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);
        w_opt=w_rand(:,i);
    end
end

gamma_l_k=(n_l(tgt_index).*abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_opt)^2)/...
    real(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);
IterVal(iter)=gamma_l_k; fprintf('Solved, Iter=%d\n',iter);

if abs(gamma_l_k-gamma_l_k_1)<e
    IterVal(iter)=gamma_l_k;
    save('mat\sinr.mat','IterVal');
    break;
end

w_k_1=w_opt; iter=iter+1;

%%%%%%NOMA-->OMA%%%%%%
% % Q_tmp
% n_l_temp=n_l;
% n_l_temp(tgt_index)=[];
% C_temp=C;
% C_temp(:,:,tgt_index)=[];
% Q_part_1_tmp=zeros(N_r*N_t,N_r*N_t);
% for l=1:L-1
%     Q_part_1_tmp=Q_part_1_tmp+n_l_temp(l)*kron(eye(N_t),C(:,:,l))'*(u_l*u_l')*kron(eye(N_t),C(:,:,l));
% end
% Q_tmp=Q_part_1_tmp+kron(eye(N_t),H_RSI)'*(u_l*u_l')*kron(eye(N_t),H_RSI)+( n*(u_l'*u_l)+P_u*sum_P_H )/P_t;
% 
% % c_m_u_tmp
% c_m_u_tmp=zeros(1,m);
% gamma_m_u_tmp=[10 8 6 4]; 
% P_u_V_H_tmp=zeros(1,M);
% for m=1:M
%     for tilde_m=1:M
%         if tilde_m~=m
%             P_u_V_H(:,m)=P_u_V_H(:,m)+P_u*abs(v_m(:,tilde_m)'*h_u(:,m))^2;
%         end
%     end
% end
% for m=1:M
%     c_m_u_tmp(:,m)=(P_u*abs(v_m(:,m)'*h_u(:,m))^2)/gamma_m_u_tmp(m)-P_u_V_H_tmp(:,m)-n_u*(v_m(:,m)'*v_m(:,m));
% end

%%%%%%Comparison Scheme%%%%%%
% cvx_begin SDP quiet
% cvx_precision high
% variable Z_tmp(N_t^2, N_t^2) hermitian
% variable y_tmp
% maximize real(trace(O*Z_tmp))
% subject to
% real(trace(Q_tmp*Z_tmp))==1;
% for k=1:K
%     real(trace(tilde_H_k(:,:,k)*Z_tmp))>=y_tmp*c_k_d(k);
% end
% for m=1:M
%     real(trace(P_m(:,:,m)*Z_tmp))<=y_tmp*c_m_u_tmp(m);
% end
% trace(Z_tmp)<=y_tmp*P_t;
% Z_tmp == semidefinite(N_t^2);
% y_tmp>=0;
% cvx_end
% 
% % gaussian randomization -- do randomization to get our feasible beamforming vectors w_cand
% tilde_W=Z_tmp/y_tmp;
% nRand = 100;
% w_rand = zeros(N_t^2, nRand);
% 
% % generate and scale to meet power constraints and other constraints
% for L_rand = 1:nRand
%     zeta(:,L_rand) = mvnrnd(zeros(N_t^2,1),tilde_W) + 1i*mvnrnd(zeros(N_t^2,1),tilde_W);         % generate nRand
%     w_rand(:,L_rand) = sqrt(trace(tilde_W))*zeta(:,L_rand)/sqrt(zeta(:,L_rand)'*zeta(:,L_rand)); % scale them so it adheres to power constraint
% end
% 
% for i=1:nRand
%     if ( ((n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
%             norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l)>maxSINR) )
%         maxSINR=(n_l(tgt_index).*norm(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
%             norm(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);
%         w_opt=w_rand(:,i);
%     end
% end
% 
% gamma_l_k=(n_l(tgt_index).*abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_opt)^2)/...
%     real(u_l'*(tilde_R_e+tilde_R_n+n.*eye(N_r*N_t))*u_l);
% IterVal(iter)=gamma_l_k; fprintf('Solved, Iter=%d\n',iter);
% 
% if abs(gamma_l_k-gamma_l_k_1)<e
%     IterVal(iter)=gamma_l_k;
%     save('mat\sinr.mat','IterVal');
%     break;
% end
% 
% w_k_1=w_opt; iter=iter+1;

end

%% Plot Figure
% % beampattern_tx
% % downlink isac
% W_opt=reshape(w_opt,N_t,N_t);
% TxBp = zeros(size(angleSpace));
% for i = 1:length(angleSpace)
%     TxBp(i) = abs(a_tx(:,i)' * (W_opt*W_opt') * a_tx(:,i))/trace(W_opt*W_opt');
% end
% figure; 
% TxBp_l = plot(angleSpaceDeg, mag2db(TxBp), 'LineWidth', 1.0, 'Linestyle', '-'); 
% for l=1:L
%     if l==tgt_index
%         tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
%     else
%         inf_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','--');
%     end
% end
% for k =1:K
%     DUE_l = line([DUE_agl(k),DUE_agl(k)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
% end
% hold on;
% grid on;
% xlabel('Angle Space [-90^\circ,90^\circ]');
% ylabel('Beampattern (dB)');
% legend([TxBp_l, tgt_l, DUE_l],'ISAC Tranmit Beampattern','Sensing Direction','Comm Direction'); 
% axis([-90, 90, min(mag2db(TxBp)), max(mag2db(TxBp))]);
% 
% % beampattern_rx
% % uplink comm
% RxBp = zeros(size(angleSpace));
% v_opt_tmp=[v_m(:,1).' v_m(:,2).' v_m(:,3).' v_m(:,4).'];
% v_opt=v_opt_tmp.';
% for i = 1:length(angleSpace)
%     RxBp(i) = abs(v_opt'*(kron(eye(M),a_rx(:,i))*kron(eye(M),a_rx(:,i))')*v_opt)/trace(v_opt*v_opt');
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.0, 'Linestyle', '-');
% for l=1:L
%     if l==tgt_index
%         tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
%     else
%         inf_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','--');
%     end
% end
% for m =1:M
%     UUE_l = line([UUE_agl(m),UUE_agl(m)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
% end
% hold on;
% grid on;
% xlabel('Angle Space [-90^\circ,90^\circ]');
% ylabel('Beampattern (dB)');
% legend([RxBp_l, tgt_l, UUE_l],'Uplink Communication Beampattern', 'Sensing Direction', 'Comm Direction'); 
% axis([-90, 90, min(mag2db(RxBp)), max(mag2db(RxBp))]);
% 
% % uplink sensing
% RxBp = zeros(size(angleSpace));
% u_opt=u_l;
% for i = 1:length(angleSpace)
%     RxBp(i) = abs(u_opt'*(kron(eye(N_t),a_rx(:,i))*kron(eye(N_t),a_rx(:,i))')*u_opt)/trace(u_opt*u_opt');
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.0, 'Linestyle', '-');
% for l=1:L
%     if l==tgt_index
%         tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
%     else
%         inf_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','--');
%     end
% end
% for m =1:M
%     UUE_l = line([UUE_agl(m),UUE_agl(m)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
% end
% hold on;
% grid on;
% xlabel('Angle Space [-90^\circ,90^\circ]');
% ylabel('Beampattern (dB)');
% legend([RxBp_l, tgt_l, UUE_l],'Sensing Beampattern', 'Sensing Direction', 'Comm Direction'); 
% axis([-90, 90, min(mag2db(RxBp)), max(mag2db(RxBp))]);
% 
% % uplink isac
% RxBp = zeros(size(angleSpace));
% v_opt_tmp=[v_m(:,1).' v_m(:,2).' v_m(:,3).' v_m(:,4).'];
% v_opt=v_opt_tmp.';
% u_opt=u_l;
% for i = 1:length(angleSpace)
%     RxBp(i) = abs(v_opt'*(kron(eye(M),a_rx(:,i))*kron(eye(M),a_rx(:,i))')*v_opt)/trace(v_opt*v_opt')+...
%         abs(u_opt'*(kron(eye(N_t),a_rx(:,i))*kron(eye(N_t),a_rx(:,i))')*u_opt)/trace(u_opt*u_opt');
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.0, 'Linestyle', '-');
% for l=1:L
%     if l==tgt_index
%         tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
%     else
%         inf_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 0.5,'linestyle','--');
%     end
% end
% for m =1:M
%     UUE_l = line([UUE_agl(m),UUE_agl(m)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
% end
% hold on;
% grid on;
% xlabel('Angle Space [-90^\circ,90^\circ]');
% ylabel('Beampattern (dB)');
% legend([RxBp_l, tgt_l, UUE_l],'Sensing Beampattern', 'Sensing Direction', 'Comm Direction'); 
% axis([-90, 90, min(mag2db(RxBp)), max(mag2db(RxBp))]);
% 
% % sum rate
% % downlink comm sinr
% gamma_d=zeros(1,K);
% h_w=zeros(1,K);
% for tilde_k=1:K
%     for k=1:K
%         if k~=tilde_k
%             h_w(tilde_k)=h_w(tilde_k)+abs(h_d(:,k)'*J_i(k)*w_opt)^2;
%         end
%     end
% end
% for k=1:K
%     gamma_d(k)=abs(h_d(:,k)'*J_i(k)*w_opt) / ( h_w(k)+P_u*norm(g(:,k))^2+n_d );
% end
% 
% % uplink comm sinr
% gamma_u=zeros(1,M);
% P_H=zeros(N_r,N_r,M);
% for tilde_m=1:M
%     if tilde_m==M
%         P_H(:,:,tilde_m)=zeros(N_r,N_r);
%     end
%     for m=tilde_m+1:M
%         P_H(:,:,tilde_m)=P_H(:,:,tilde_m)+P_u*( h_u(:,m)*h_u(:,m)' );
%     end
% end
% for m=1:M
%     gamma_u(m)=P_u*abs(v_m(:,m)'*h_u(:,m))^2 / real( v_m(:,m)' * (P_H(:,:,m)+R_e+R_n+n_u*eye(N_r)) * v_m(:,m) );
% end