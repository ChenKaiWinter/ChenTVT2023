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

% transmit antennas steering vecter all angl
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
load('mat\a_tx.mat');
load('mat\a_rx.mat');
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

%% First Slot
P_u_H_u=zeros(N_r,N_r,M);
v_m=zeros(N_r,M);
for tilde_m=1:M
    if tilde_m==M
        P_u_H_u(:,:,tilde_m)=zeros(N_r,N_r);
    end
    for m=tilde_m+1:M
        P_u_H_u(:,:,tilde_m)=P_u_H_u(:,:,tilde_m)+P_u*(h_u(:,m)*h_u(:,m)');
    end
end
for m=1:M
    v_m(:,m)=sqrt(P_u)*(P_u_H_u(:,:,m)+n_u*eye(N_r))^(-1)*h_u(:,m);
end

%% Second Slot
while true

% Iter I - u_l: the sensing beamforming vector
u_l=( (tilde_R_e+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1 )...
    /(w_k_1'*kron(eye(N_t),C(:,:,tgt_index))'*(tilde_R_e+n*eye(N_r*N_t))^(-1)*kron(eye(N_t),C(:,:,tgt_index))*w_k_1);

gamma_l_k_1=(n_l(tgt_index)*(abs(u_l'*kron(eye(N_t),C(:,:,tgt_index))*w_k_1)^2))/...
    real(u_l'*(tilde_R_e+n.*eye(N_r*N_t))*u_l);

% Iter II - w: the transmit beamforming vector
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
Q=Q_part_1+n*(u_l'*u_l)/P_t;

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
    c_k_d(:,k)=n_d;
end

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
trace(Z)<=y*P_t;
Z == semidefinite(N_t^2);
y>=0;
cvx_end

% eigenvalue decomposition: sqrt(dominant eigenvalue)*dominant eigenvector (It has been proved that there is no rank one solution)
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

end

%% Plot Figure
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
%     gamma_d(k)=abs(h_d(:,k)'*J_i(k)*w_opt) / ( h_w(k)+n_d );
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
%     gamma_u(m)=P_u*abs(v_m(:,m)'*h_u(:,m))^2 / real( v_m(:,m)' * (P_H(:,:,m)+n_u*eye(N_r)) * v_m(:,m) );
% end