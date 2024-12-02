function data = gmac_EbN0(k, Lmin, Lmax, alpha, ...
    EbN0db_lower, EbN0db_upper, adjustRadii, fixP1)
% This script computes and plots active user spectral efficiency v EbN0_dB. 
% Active user spectral efficiency=k*E[Ka]/n_total = L(1-alpha)/n via L, 
% keeping alpha and n fixed. Note unlike in the paper, alpha in the code
% represents the probability each user is silent, not active. 
% 
% Specifically, for the number of users L in [Lmin, Lmax], find the minimum 
% EbN0_dB needed in [EbN0db_lower, EbN0db_upper] to achieve a specified 
% error, fixing alpha.
% 
% Setting adjustRadii=true uses the minimum r that can attain
% target_epsTotal.
% Setting fixP1=true skips the golden section search over P1.

addpath RCU_KaUnknown_SRA;


myPool = parpool(8); % num workers
fprintf('Running gmac_EbN0...\n')
tStart = tic;
n = 2000; % total num real channel uses
num_L = 20;
L_arr = ceil(linspace(Lmin, Lmax, num_L));
var_arr = L_arr*alpha*(1-alpha); 
rad_l_arr = ceil(var_arr); 
rad_u_arr = rad_l_arr;  
target_epsTotal = 1e-1;
eff_mu_arr = L_arr*(1-alpha) / n; % effective user density: E[Ka]/n
bin_search_max_iter = 40; % for binary search

% When the target_epsTotal is achieved, epsMD,epsFA,epsAUE and the minimum 
% Eb/N0 needed are:
epsMD = zeros(num_L,1);
epsFA = zeros(num_L,1);
epsAUE = zeros(num_L,1); 
epsTotal = zeros(num_L,1);
optP1     = zeros(num_L,1);
min_EbN0db = zeros(num_L,1);
bin_search_num_iter_conv = zeros(num_L,1); % up till binary search converges

% Check whether epsTotal decreases with P:
bin_search_epsTotal_decreases_w_P = false(num_L,1);
% Investigate relationship between optP1 and P via binary search history:
bin_search_P = zeros(num_L,bin_search_max_iter+2);
bin_search_optP1 = zeros(num_L,bin_search_max_iter+2);

for iL = 1:num_L
    L = L_arr(iL);
    rad_l = rad_l_arr(iL);
    rad_u = rad_u_arr(iL);
    fprintf('L=%d [%d/%d]\n', L, iL, num_L);
    data0 = EbN0_KaBinomialUnknown_SRA(k, n, L, alpha, target_epsTotal, ...
        rad_l, rad_u, EbN0db_lower, EbN0db_upper, ...
        adjustRadii, fixP1);
    epsMD(iL) = data0.epsMD_RCU;
    epsFA(iL) = data0.epsFA_RCU;
    epsAUE(iL) = data0.epsAUE_RCU;
    epsTotal(iL) = data0.epsTotal_RCU;
    optP1(iL) = data0.optP1;
    min_EbN0db(iL) = data0.EbN0db_RCU;
    % Binary search data:
    num_iter_conv = data0.bin_search_data.num_iter_conv;
    bin_search_num_iter_conv(iL) = num_iter_conv;
    bin_search_epsTotal_decreases_w_P(iL) = data0.bin_search_data.epsTotal_decreases_w_P;
    bin_search_P(iL, 1:num_iter_conv+2) = data0.bin_search_data.P_arr;
    bin_search_optP1(iL, 1:num_iter_conv+2) = data0.bin_search_data.optP1_arr;
end
fprintf('[Reached the end of gmac_EbN0 in %.2f]', toc(tStart));


data.k = k;
data.n = n;
data.L = L_arr;
data.alpha = alpha;
data.p_Ka = 'Binomial';
data.rad_l_arr = rad_l_arr;
data.rad_u_arr = rad_u_arr;
data.adjustRadii = adjustRadii;
data.fixP1 = fixP1;
data.target_epsTotal = target_epsTotal;
data.obj = 'max_pTotal';

data.eff_mu = eff_mu_arr;
data.epsMD = epsMD;
data.epsFA = epsFA;
data.epsAUE = epsAUE;
data.optP1 = optP1;
data.min_EbN0db = min_EbN0db;

data.bin_search_num_iter_conv = bin_search_num_iter_conv;
data.bin_search_epsTotal_decreases_w_P = bin_search_epsTotal_decreases_w_P;
data.bin_search_P = bin_search_P;
data.bin_search_optP1 = bin_search_optP1;

alpha_str = sprintf('%.1f,', data.alpha);
paramStr = ['obj=' data.obj '_' ...
    sprintf('target_epsTotal=%0.1fx10^%i', ...
    10^mod(log10(data.target_epsTotal),1), ...
    floor(log10(data.target_epsTotal))) ...
    sprintf('k=%d,n=%d,alpha=%s', data.k, data.n, alpha_str)];
dt = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss');
dtStr = char(dt);
filename = ['mu_EbN0_' paramStr dtStr];
save([filename '.mat'], 'data', '-v7.3');

fprintf(['\nIn binary search: epsTotal doesnt decrease monotonically with P ' ...
    'in %d/%d searches'], ...
    sum(~bin_search_epsTotal_decreases_w_P), num_L);

figure;
plot(data.min_EbN0db, data.eff_mu, 'marker', '.', 'LineWidth', 2);
hold on;
xlabel('Eb/N0 (dB)');
ylabel('Effective user density mu=E[Ka]/n');
fontSize = 18;
set(gca, 'FontSize', fontSize);
saveas(gcf, [filename '.fig']);  
hold off;
delete(myPool);
end