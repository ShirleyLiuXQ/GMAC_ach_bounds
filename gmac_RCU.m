function data = gmac_RCU(P1_asFactorOfP)
% This script computes and plots pErrors v EbN0_dB
addpath RCU_KaUnknown_SRA;

save_path='results/';
dt = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss');
dtStr = char(dt);
log_file_name = [save_path 'pError_v_EbN0_' dtStr '.txt'];

fprintf('Running gmac_RCU...\n')
tStart = tic;
k = 8; % Each codeword carries k bits. Each codebook contains M=2^k codewords.
% n = 2000; % total num real channel uses = amp_n * k

EbN0db_lower = 0;   
EbN0db_upper = 20;  
num_EbN0 = 10;
EbN0db_vec = linspace(EbN0db_lower, EbN0db_upper, num_EbN0);
amp_n_over_L = 5;
L = 50;  
amp_n = round(amp_n_over_L * L);
n = amp_n * k; % total number of channel uses
alpha = [0.1, 0.5, 0.9]; % 0.5; % alpha denotes fraction of silent users
rad_l = 0; %[0, 1, Inf];
rad_u = 2; %[0, 1, Inf];
tail_prob=1e-13;
frac_TOL_golden=0.01; % TOL for golden section search over P1\in[0,P]. 
% Finer grids give smoother pErr vs EbN0 plots

fprintf(['Parameters: k=%d, n=%d, L=%d, alpha=%.2f, EbN0db_lower=%.2f, ' ...
    'EbN0db_upper=%.4f, rad_l=%.2f, rad_u=%.2f, tail_prob=%s, frac_TOL_golden=%.2f\n'], ...
    k, n, L, alpha, EbN0db_lower, EbN0db_upper, ...
    rad_l, rad_u, sprintf('%.0e', tail_prob), frac_TOL_golden);
if isempty(P1_asFactorOfP)
    obj = 'max_pTotal'; % 'max'
else
    obj = [];
end
data = RCU_KaBinomialUnknown_SRA(k, n, L, alpha, EbN0db_vec, rad_l, rad_u, ...
    tail_prob, obj, P1_asFactorOfP, frac_TOL_golden, log_file_name);
fprintf('[Reached the end of std_gmac_RCU in %.2f]', toc(tStart));

rad_l_str = sprintf('%d,', data.rad_lower);
rad_u_str = sprintf('%d,', data.rad_upper);
alpha_str = sprintf('%.1f,', data.alpha);

if isempty(P1_asFactorOfP)
    paramStr = ['obj=' data.obj, '_' sprintf('k=%d,n=%d,L=%d,alpha=%srl=%sru=%stail=%s,P1TOL=%.2f_', ...
        data.k, data.n, data.L, alpha_str, rad_l_str, rad_u_str, sprintf('%.0e', tail_prob), frac_TOL_golden)];
else
    paramStr = [sprintf('P1factor=%.2f',P1_asFactorOfP), ...
        '_' sprintf('k=%d,n=%d,L=%d,alpha=%srl=%sru=%stail=%s,P1TOL=%.3f_', ...
        data.k, data.n, data.L, alpha_str, rad_l_str, rad_u_str, sprintf('%.0e', tail_prob), frac_TOL_golden)];
end

dt = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss');
dtStr = char(dt);
filename = ['pErr_EbN0_' paramStr dtStr];
save([filename '.mat'], 'data', '-v7.3');

%% 
figure;
blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
plotMD = plot(data.EbN0db, data.pMD(:,1), 'marker', 'o', 'LineWidth', 2, 'Color', blue); 
hold on;
plotFA = plot(data.EbN0db, data.pFA(:,1), 'marker', '+', 'LineWidth', 2, 'Color', red);
plotAUE = plot(data.EbN0db, data.pAUE(:,1), 'marker', '*', 'LineWidth', 2, 'Color', yellow);

plot(data.EbN0db, data.pMD(:,2), 'marker', 'o', 'LineWidth', 2, 'Color', blue, 'LineStyle','--'); 
plot(data.EbN0db, data.pFA(:,2), 'marker', '+', 'LineWidth', 2, 'Color', red, 'LineStyle','--');
plot(data.EbN0db, data.pAUE(:,2), 'marker', '*', 'LineWidth', 2, 'Color', yellow, 'LineStyle','--');   
plot(data.EbN0db, data.pMD(:,3), 'marker', 'o', 'LineWidth', 2, 'Color', blue, 'LineStyle',':'); 
plot(data.EbN0db, data.pFA(:,3), 'marker', '+', 'LineWidth', 2, 'Color', red, 'LineStyle',':');
plot(data.EbN0db, data.pAUE(:,3), 'marker', '*', 'LineWidth', 2, 'Color', yellow, 'LineStyle',':');   
set(gca, 'YScale', 'log');

plot(data.EbN0db, data.floor_pMD(1) * ones(num_EbN0, 1), 'Color', blue)
plot(data.EbN0db, data.floor_pFA(1) * ones(num_EbN0, 1), 'Color', red)
plot(data.EbN0db, data.floor_pAUE(1) * ones(num_EbN0, 1), 'Color', yellow)
plot(data.EbN0db, data.floor_pMD(2) * ones(num_EbN0, 1), 'Color', blue, 'LineStyle', '--')
plot(data.EbN0db, data.floor_pFA(2) * ones(num_EbN0, 1), 'Color', red, 'LineStyle', '--')
plot(data.EbN0db, data.floor_pAUE(2) * ones(num_EbN0, 1), 'Color', yellow, 'LineStyle', '--')
plot(data.EbN0db, data.floor_pMD(3) * ones(num_EbN0, 1), 'Color', blue, 'LineStyle', ':')
plot(data.EbN0db, data.floor_pFA(3) * ones(num_EbN0, 1), 'Color', red, 'LineStyle', ':')
plot(data.EbN0db, data.floor_pAUE(3) * ones(num_EbN0, 1), 'Color', yellow, 'LineStyle', ':')

plot_r1 = plot(NaN,NaN, 'Color', 'k', 'LineWidth', 2);
plot_r2 = plot(NaN,NaN, 'Color', 'k', 'LineWidth', 2, 'LineStyle','--');
plot_r3 = plot(NaN,NaN, 'Color', 'k', 'LineWidth', 2, 'LineStyle',':');

if isscalar(data.rad_lower) && isscalar(data.rad_upper)
    % different alpha, fixed r
    legend([plotMD, plotFA, plotAUE, plot_r1, plot_r2, plot_r3], ...
        'p_M_D', 'p_F_A', 'p_A_U_E', ...
        sprintf('alpha=%.1f', data.alpha(1)), ...
        sprintf('alpha=%.1f', data.alpha(2)), ...
        sprintf('alpha=%.1f', data.alpha(3)), 'Location', 'southwest');
elseif isscalar(data.alpha)
    % different r, fixed alpha
    legend([plotMD, plotFA, plotAUE, plot_r1, plot_r2, plot_r3], ...
        'p_M_D', 'p_F_A', 'p_A_U_E', ...
        sprintf('r=%d', data.rad_lower(1)), ...
        sprintf('r=%d', data.rad_lower(2)), ...
        sprintf('r=%d', data.rad_lower(3)));
end

% title(param_str);
xlabel('Eb/N0 (dB)');
ylabel('MD, FA, AUE probabilities');
fontSize = 18;
set(gca, 'FontSize', fontSize);
figname = ['pErr_EbN0_' paramStr dtStr];
saveas(gcf, [figname '.fig']);  
hold off;

end