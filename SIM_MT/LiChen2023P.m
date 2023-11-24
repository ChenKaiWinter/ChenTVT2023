% "Signal Processing and Beamforming Optimization for Full-Duplex NOMA Enabled Integrated Sensing and Communication" 
% - Kai Chen et. al.
% Kai Chen, 10/07/2023

%% Set up the System -- Multi Target
K = 2; DUE_dis = [128 128]; DUE_agl = [-50 50];           % K: number of DUE
M = 4; UUE_dis = [8 16 32 64]; UUE_agl = [-70 -30 40 60]; % M: number of UUE
L = 3; tgt_agl = [-20 0 20];                              % L: number of tgt

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

%% Algorithm
% w_intial
w_0 = sqrt(P_t/(2*N_t^2))*randn(N_t^2, 1)+sqrt(P_t/(2*N_t^2))*1i*randn(N_t^2, 1);
save('mat\w_0.mat','w_0');
load('mat\w_0.mat'); % get the inital power

% ch_intial
G_dB = 10;                                        % Rician factor
beta = 3.8;
[a_rx, a_tx]=sv(N_t,N_r,lambda,angleSpace,TxAntLoc,RxAntLoc,DUE_agl,UUE_agl,tgt_agl);
[h_d, h_u, g, H_RSI] = gen_h_rci(G_dB, beta, a_tx, a_rx, N_t, N_r, K, M, ...
    DUE_agl, UUE_agl, DUE_dis, UUE_dis, D_U_dis); % h_d, h_u, g are chanel matrix
save('mat\ch.mat','h_d','h_u','g','H_RSI');       % save the channel
load('mat\ch.mat');                               % get the inital chanel

% generate the matrix that I need
w_k_1=w_0; % w^(k-1)
[C, J_i, R_e, R_n, tilde_R_e, tilde_R_n, R_k, sum_R_i] = gen_mat(N_r, N_t, K, L, H_RSI, a, b, n_l, w_k_1);

% begin
e = 0.001; iter = 1; maxSINR = 0;

%% Iterative Algorithm
% IBS inital
zeta_l=0;
zeta_u=1e2;
zeta=(zeta_l+zeta_u)/2;

while true

% Iter I - u_l: the sensing beamforming vector
u_l=zeros(N_r*N_t,L);
for l=1:L
    u_l(:,l)=( (tilde_R_e(:,:,l)+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,l))*w_k_1 )...
    /(w_k_1'*kron(eye(N_t),C(:,:,l))'*(tilde_R_e(:,:,l)+tilde_R_n+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,l))*w_k_1);
end

gamma_l_k_1=zeros(1,L);
for l=1:L
    gamma_l_k_1(l)=(n_l(l)*(abs(u_l(:,l)'*kron(eye(N_t),C(:,:,l))*w_k_1)^2))/...
    norm(u_l(:,l)'*(tilde_R_e(:,:,l)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,l));
end

% Iter II - v_m: the uplink receive communication beamforming vector
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
    v_m(:,m)=sqrt(P_u)*( (P_u_H_u(:,:,m)+R_e+R_n+n_u*eye(N_r))^(-1) )*h_u(:,m);
end

% Iter III - w: the transmit beamforming vector
% O_l
O=zeros(N_t^2,N_t^2,L);
for l=1:l
    O(:,:,l)=n_l(l)*(kron(eye(N_t),C(:,:,l))'*(u_l(:,l)*u_l(:,l)')*kron(eye(N_t),C(:,:,l)));
end

% tilde_Q_l
tilde_Q=zeros(N_t^2,N_t^2,L);
for tgt_index=1:L
    n_l_temp=n_l;
    n_l_temp(tgt_index)=[];
    C_temp=C;
    C_temp(:,:,tgt_index)=[];
    Q_part_1=zeros(N_r*N_t,N_r*N_t);
    for l=1:L-1
        Q_part_1=Q_part_1+n_l_temp(l)*kron(eye(N_t),C(:,:,l))'*(u_l(:,tgt_index)*u_l(:,tgt_index)')*kron(eye(N_t),C(:,:,l));
    end
    tilde_Q(:,:,tgt_index)=Q_part_1+kron(eye(N_t),H_RSI)'*(u_l(:,tgt_index)*u_l(:,tgt_index)')*kron(eye(N_t),H_RSI)...
        +n.*(u_l(:,tgt_index)'*u_l(:,tgt_index))/P_t;
end

% tilde_H_k
gamma_k_d=[18 18];
tilde_H_k_part_2=zeros(N_t^2,N_t^2,K);
for k=1:K
    h_d_temp=h_d;
    h_d_temp(:,k)=[];
    for k_temp=1:K-1
        tilde_H_k_part_2(:,:,k)=tilde_H_k_part_2(:,:,k)+J_i(k_temp)'*h_d_temp(:,k_temp)*h_d_temp(:,k_temp)'*J_i(k_temp);
    end
end
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
            P_m_part_2(:,:,M)=n_l(l).*(J_i(i)'*C(:,:,l)'*v_m(:,m)*v_m(:,m)'*C(:,:,l)*J_i(i));
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

% use CVX to solve Iter III
cvx_begin SDP quiet
cvx_precision high
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%variable%%%
variable tilde_W(N_t^2, N_t^2) hermitian
variable y
%%%%OBJ_FX%%%%
% This is a problem with no objective function
%%constraint%%
subject to
%%%%%%constraint(59b)%%%%%%
for l=1:L
    real(trace(O(:,:,l)*tilde_W))-zeta*(real(trace(tilde_Q(:,:,l)*tilde_W)))>=0;
end
%%%%%%constraint(59c)%%%%%%
for k=1:K
    real(trace(tilde_H_k(:,:,k)*tilde_W))-y*c_k_d(k)>=0;
end
%%%%%%constraint(59d)%%%%%%
for m=1:M
    y*c_m_u(m)-real(trace(P_m(:,:,m)*tilde_W))>=0;
end
%%%%%%constraint(59e)%%%%%%
trace(tilde_W)<=P_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tilde_W==semidefinite(N_t^2);
y>=0;
cvx_end

%%%%%%solution%%%%%%
% Gaussian randomization -- do randomization to get our feasible beamforming vectors w_cand
nRand = 100;
w_rand = zeros(N_t^2, nRand);

% generate and scale to meet power constraints and other constraints
for L_rand = 1:nRand
    % generate nRand
    tmp(:,L_rand) = mvnrnd(zeros(N_t^2,1),tilde_W) + 1i*mvnrnd(zeros(N_t^2,1),tilde_W);
    % scale them so it adheres to power constraint
    w_rand(:,L_rand) = sqrt(trace(tilde_W))*tmp(:,L_rand)/sqrt(tmp(:,L_rand)'*tmp(:,L_rand));
end

for i=1:nRand
    % find the min SINR and corresponding vector
    minSINR=(n_l(1).*norm(u_l(:,1)'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm( u_l(:,1)'*(tilde_R_e(:,:,1)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,1) );
    minu_l=u_l(:,1);
    if( ((n_l(2).*norm(u_l(:,2)'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l(:,2)'*(tilde_R_e(:,:,2)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,2)))<minSINR )
        minSINR=((n_l(2).*norm(u_l(:,2)'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l(:,2)'*(tilde_R_e(:,:,2)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,2)));
        minu_l=u_l(:,2);
    end
    if( ((n_l(3).*norm(u_l(:,3)'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l(:,3)'*(tilde_R_e(:,:,3)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,3)))<minSINR )
        minSINR=((n_l(2).*norm(u_l(:,2)'*kron(eye(N_t),C(:,:,tgt_index))*w_rand(:,i))^2)/...
            norm(u_l(:,2)'*(tilde_R_e(:,:,3)+tilde_R_n+n.*eye(N_r*N_t))*u_l(:,3)));
        minu_l=u_l(:,2);
    end
    if minSINR>maxSINR 
        maxSINR=minSINR;
        w_opt=w_rand(:,i);
    end
end

if ( (strcmp(cvx_status,'Solved')==1)||(strcmp(cvx_status,'Inaccurate/Solved')==1) ) % if there is a solution
    zeta_l=zeta;
    zeta=(zeta_l+zeta_u)/2;
    fprintf('Solved, Upward, zeta=%.5f, zeta_l=zeta, (Iter: zeta_l=%.5f, zeta_u=%.5f, zeta_u-zeta_l=%.5f)\n',zeta,zeta_l,zeta_u,zeta_u-zeta_l);
    if (zeta_u-zeta_l<0.1)
        break;
    end
else
    zeta_u=zeta;
    zeta=(zeta_l+zeta_u)/2;
    fprintf('Unresolved, Downward, zeta=%.5f, zeta_u=zeta, (Iter: zeta_l=%.5f, zeta_u=%.5f, zeta_u-zeta_l=%.5f)\n',zeta,zeta_l,zeta_u,zeta_u-zeta_l);
end

zetaIter(iter)=zeta;
save('mat\SensingSINR.mat','zetaIter');

w_k_1=w_opt;iter=iter+1;

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
% TxBp_l = plot(angleSpaceDeg, mag2db(TxBp), 'LineWidth', 1.5, 'Linestyle', '-'); 
% for l=1:L
%     tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(TxBp)), max(mag2db(TxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
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
% v_opt=120*v_opt_tmp.';
% for i = 1:length(angleSpace)
%      RxBp(i)=abs(v_opt'*(kron(eye(M),a_rx(:,i))*kron(eye(M),a_rx(:,i))')*v_opt)/trace(v_opt*v_opt');
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.5, 'Linestyle', '-');
% for l=1:L
%     tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
% end
% for m =1:M
%     UUE_l = line([UUE_agl(m),UUE_agl(m)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'magenta', 'LineWidth', 1.0,'linestyle','--');
% end
% hold on;
% grid on;
% xlabel('Angle Space [-90^\circ,90^\circ]');
% ylabel('Beampattern (dB)');
% legend([RxBp_l, tgt_l, UUE_l],'Uplink Receive Beampattern', 'Sensing Direction', 'Comm Direction'); 
% axis([-90, 90, min(mag2db(RxBp)), max(mag2db(RxBp))]);
% 
% % uplink sensing
% RxBp = zeros(size(angleSpace));
% u_opt=u_l;
% for i = 1:length(angleSpace)
%     RxBp(i) = abs(u_opt'*(kron(eye(N_t),a_rx(:,i))*kron(eye(N_t),a_rx(:,i))')*u_opt)//trace(u_opt*u_opt')
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.5, 'Linestyle', '-');
% for l=1:L
%     tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
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
% u_opt=zeros(N_r*N_t,1);
% for l=1:L
%     u_opt=u_opt+u_l(:,l);
% end
% for i = 1:length(angleSpace)
%     RxBp(i) = abs(v_opt'*(kron(eye(M),a_rx(:,i))*kron(eye(M),a_rx(:,i))')*v_opt)/trace(v_opt*v_opt')+...
%         abs(u_opt'*(kron(eye(N_t),a_rx(:,i))*kron(eye(N_t),a_rx(:,i))')*u_opt)/trace(u_opt*u_opt');
% end
% figure; 
% RxBp_l = plot(angleSpaceDeg, mag2db(RxBp), 'LineWidth', 1.5, 'Linestyle', '-');
% for l=1:L
%     tgt_l = line([tgt_agl(l),tgt_agl(l)],[min(mag2db(RxBp)), max(mag2db(RxBp))],'Color', 'black', 'LineWidth', 1.0,'linestyle','--');
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