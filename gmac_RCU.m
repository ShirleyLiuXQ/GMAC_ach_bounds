function data = gmac_RCU(P1_asFactorOfP)
% This script computes and plots pErrors v EbN0_dB.

addpath RCU_KaUnknown_SRA;

fprintf('Running gmac_RCU...\n')
tStart = tic;
k = 6; % Each codeword carries k bits. Each codebook contains M=2^k codewords.


EbN0db_lower = 10; 
EbN0db_upper = 15;
num_EbN0 = 10;
EbN0db_vec = linspace(EbN0db_lower, EbN0db_upper, num_EbN0);
amp_n_over_L = 0.45;
L = 500; 
amp_n = round(amp_n_over_L * L);
n = amp_n * k; % total number of channel uses
alpha = 0.3; % alpha is the probability that a user is inactive (opposite
% to the convention in our paper)
rad_l = 10;  
rad_u = 10; 
if isempty(P1_asFactorOfP)
    obj = 'max_pTotal'; % 'max'
else
    obj = [];
end

data = RCU_KaBinomialUnknown_SRA(k, n, L, alpha, EbN0db_vec, rad_l, rad_u, obj, P1_asFactorOfP);
fprintf('[Reached the end of gmac_RCU in %.2f]', toc(tStart));

rad_l_str = sprintf('%d,', data.rad_lower);
rad_u_str = sprintf('%d,', data.rad_upper);
alpha_str = sprintf('%.1f,', data.alpha);
if isempty(P1_asFactorOfP)
    paramStr = ['obj=' data.obj, '_' sprintf('k=%d,n=%d,L=%d,alpha=%srl=%sru=%s', ...
        data.k, data.n, data.L, alpha_str, rad_l_str, rad_u_str)];
else
    paramStr = [sprintf('P1factor=%.2f',P1_asFactorOfP), ...
        '_' sprintf('k=%d,n=%d,L=%d,alpha=%srl=%sru=%s', ...
        data.k, data.n, data.L, alpha_str, rad_l_str, rad_u_str)];
end
dt = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss');
dtStr = char(dt);
filename = ['pErr_EbN0_' paramStr dtStr];
save([filename '.mat'], 'data', '-v7.3');

figure;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
plotMD = plot(data.EbN0db, data.pMD(:,1), 'marker', 'o', 'LineWidth', 2, 'Color', blue); 
hold on;
plotFA = plot(data.EbN0db, data.pFA(:,1), 'marker', '+', 'LineWidth', 2, 'Color', red);
plotAUE = plot(data.EbN0db, data.pAUE(:,1), 'marker', '*', 'LineWidth', 2, 'Color', yellow); 
set(gca, 'YScale', 'log');

plot(data.EbN0db, data.floor_pMD(1) * ones(num_EbN0, 1), 'Color', blue)
plot(data.EbN0db, data.floor_pFA(1) * ones(num_EbN0, 1), 'Color', red)
plot(data.EbN0db, data.floor_pAUE(1) * ones(num_EbN0, 1), 'Color', yellow)

xlabel('Eb/N0 (dB)');
ylabel('MD, FA, AUE probabilities');
fontSize = 18;
set(gca, 'FontSize', fontSize);
figname = ['pErr_EbN0_' paramStr dtStr];
saveas(gcf, [figname '.fig']);  
if ~isRemote
    export_fig(gcf, figname, '-nocrop', '-pdf', '-m2', '-transparent', '-q101'); 
end
hold off;

end