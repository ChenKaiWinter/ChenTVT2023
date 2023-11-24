K = 2; DUE_dis = [128 128]; DUE_agl = [-50 50];           % K: number of DUE
M = 4; UUE_dis = [8 16 32 64]; UUE_agl = [-70 -30 40 60]; % D: number of UUE
L = 3; tgt_agl = [-10 0 20];                              % L: number of tgt
tgt_index=2;                                              % index of the target
for k=1:K
    polarscatter(pi/2-DUE_agl(k)/180*pi,DUE_dis(k),80,'r','filled');hold on;
end
for m=1:M
    polarscatter(pi/2-UUE_agl(m)/180*pi,UUE_dis(m),80,'b','s','filled');hold on;
end
for l=1:L
    if l==tgt_index
        polarplot([pi/2-tgt_agl(tgt_index)/180*pi,pi/2-tgt_agl(tgt_index)/180*pi],[150,0],'black--','LineWidth',1);hold on;
    else
        polarplot([pi/2-tgt_agl(l)/180*pi,pi/2-tgt_agl(l)/180*pi],[150,0],'black--');hold on;
    end
end
for k=1:K
    polarplot([pi/2-DUE_agl(k)/180*pi,pi/2-DUE_agl(k)/180*pi],[0,150],'r--','LineWidth',1);hold on;
end
for m=1:M
    polarplot([pi/2-UUE_agl(m)/180*pi,pi/2-UUE_agl(m)/180*pi],[0,150],'b--','LineWidth',1);hold on;
end
thetalim([0 180])