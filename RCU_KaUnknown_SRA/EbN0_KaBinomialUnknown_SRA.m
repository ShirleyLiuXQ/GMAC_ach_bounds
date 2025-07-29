function data = EbN0_KaBinomialUnknown_SRA(k, n, L, alpha, ...
    target_epsTotal, rad_l, rad_u, tail_prob, ...
    EbN0db_lower, EbN0db_upper, P1_asFactorOfP, ...
    frac_TOL_binary, frac_TOL_golden, adjustRadii, log_file_name)
% Find the minimal required EbN0 (in dB) such that the total error 
% probability max(pMD,pFA)+pAUE is below the prescribed threshold
% target_epsTotal, for a system whose number of active users is 
% unknown and follows a Binomial distribution.
% 
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of real DoFs)
%   L   : total num users, among which Ka is active
%   alpha : probability that a user is silent, in replacement of E_Ka
%   target_epsTotal : target total probability
%   rad_l : lower decoding radius
%   rad_u : upper decoding radius
%   The range of power to search over: [EbN0db_lower, EbN0db_upper] 
%
% OUTPUT
%   data : store the system parameters and the minimal required EbN0 in dB
% 
% This function requires target_epsTotal as input and returns the minimal 
% EbN0 needed as output, while its counterpart RCU_KaBinomialUnknown_SRA 
% requires EbN0 as input and returns the eps_MD, eps_FA, eps_AUE for that
% given EbN0.

log_file = fopen(log_file_name, 'a');
tStart = tic;
fprintf(log_file, 'Running EbN0_KaBinomialUnknown_SRA ...\n');
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k = 128; % Number of bits
    n = 19200; % number real channel uses
    L = 100;
    alpha = 0.5;
    target_epsTotal = 0.1; 
    rad_l = 0;
    rad_u = 0;
    tail_prob=1e-7;
    EbN0db_lower=8;
    EbN0db_upper=14;
    P1_asFactorOfP=0.7;
    frac_TOL_binary = 1/20;
    frac_TOL_golden = 1/20;
    adjustRadii=1;
end

%% Ka is Binomial. Can be modified into other distributions.
p_Ka = @(K) binopdf(K,L,1-alpha);
E_Ka = L*(1-alpha);

% the factor of 2 below is because N0=2sigma^2=2
P_lower = 2*k*10^(EbN0db_lower/10)/n;
P_upper = 2*k*10^(EbN0db_upper/10)/n;

% adjustRadii = 0;
% adjustRadii: set to false by default. The range of search is defined by
% rad_l,rad_u and tail_prob. When the range is small, adjustRadii=true
% increases the range to ensure the asymptotic error floors lie below the
% target total error target_epsTotal. When the range is big,
% adjustRadii=true shrinks the range to speed up runtime. However the
% adjustments can be misleading for different target_epsTotal so by default
% set to false.
if adjustRadii == 1
    % Find suitable decoding radii:
    [rad_l, rad_u] = radii_for_target_eps_SRA(n,L,tail_prob, E_Ka,p_Ka, ...
        target_epsTotal);
end
%% Function to compute the RCU bound on three errors
f_rcu = @(P,P1) RCU_KaRandomUnknown_SRA(P,P1,rad_l,rad_u, tail_prob, ...
    k,n,L,E_Ka,p_Ka, log_file_name); % function handle

%% Search for the minimal required EbN0 that gives 
% max(eps_MD,eps_FA)+eps_AUE below the specified threshold target_epsTotal.
[epsMD_RCU, epsFA_RCU, epsAUE_RCU, epsTotal_RCU, P_RCU, optP1, bin_search_data] = ...
    binary_search_P_SRA(f_rcu, P_lower, P_upper, target_epsTotal, ...
    target_epsTotal*frac_TOL_binary, P1_asFactorOfP, frac_TOL_golden, log_file_name); 
% use 1/20 TOL instead of 1/100 TOL for speed
if P_RCU ~= -1
    min_EbN0db_RCU = 10*log10(n*P_RCU/k/2);
    fprintf(log_file, ['binary search over P has converged to min_EbN0db_RCU=%e ' ...
        'epsTotal_RCU=%e\n'], min_EbN0db_RCU, epsTotal_RCU);
else
    min_EbN0db_RCU = -1; % the range [EbN0db_lower, EbN0db_upper] is too 
    % low to reach target_epsTotal
    fprintf(log_file, 'binary search over P failed to find EbN0 in range that achieves target eps_total');
end

%% Save the results
sim_time = toc(tStart);
data.k = k;
data.n = n;
data.L = L;
data.alpha = alpha;
data.target_epsTotal = target_epsTotal;
data.rad_l = rad_l;
data.rad_u = rad_u;
data.tail_prob = tail_prob;
data.adjustRadii = adjustRadii;
data.EbN0db_lower = EbN0db_lower;
data.EbN0db_upper = EbN0db_upper;
data.P1_asFactorOfP = P1_asFactorOfP;
data.frac_TOL_binary = frac_TOL_binary;
data.frac_TOL_golden = frac_TOL_golden;
data.E_Ka = E_Ka;
data.p_Ka = 'Binomial';
data.obj = 'max_pTotal';

data.epsMD_RCU = epsMD_RCU;
data.epsFA_RCU = epsFA_RCU;
data.epsAUE_RCU = epsAUE_RCU;
data.epsTotal_RCU = epsTotal_RCU;
data.optP1     = optP1;
data.bin_search_data = bin_search_data;
data.min_EbN0db_RCU = min_EbN0db_RCU;
data.sim_time = sim_time;
fprintf(log_file, '[EbN0_KaBinomialUnknown_SRA finished in %f secs]\n', sim_time);

% if DEBUG ~= 1
%     filename = ['EbN0_KaPoissonUnknown_EKa_' num2str(E_Ka) '_target_epsMD_' ...
%         num2str(target_eps_MD) '_target_epsFA_' num2str(target_eps_FA) ...
%         '_k_' num2str(k) '_n_' num2str(n) '.mat'];
%     save(filename, 'data', '-v7.3');
% else
%     keyboard
% end

end

function [rad_l, rad_u] = radii_for_target_eps_SRA(n,L, tail_prob, E_Ka,p_Ka, ...
    target_epsTotal)
% Find the suitable decoding radii rad_l,rad_u to ensure the asymptotic
% error floors lie below the target total error target_epsTotal, and that  
% one cannot significantly further reduce the minimal Eb/N0 needed to 
% achieve the target error by increasing the decoding radii, as explained 
% on p.9 [Ngo22].
% 
% The rationale behind this function is that since 
% 1) The asymptotic actual errors eps_MD,eps_FA,eps_AUE are lower bounded 
% by the asymptotic error floors.
% 2) The non-asymptotic/ finite blocklength errors are worse than the
% asymptotic errors, and thus the asymptotic error floors.
% If the asymptotic total error floor (given rad_l,rad_u) i.e.
% max(floor_epsMD,floor_epsFA)+floor_epsAUE > target_epsTotal, then the
% decoder defined with these rad_l,rad_u will never give a total error
% below the target_epsTotal.

fprintf(log_file, "Start adjusting rad_l and rad_u using tail_prob and target_epsTotal\n");
% Initialise:
rad_l = 0; 
rad_u = 0;
% Zero decoding radii give the worst decoding performance in terms of the
% metric m chosen in eq(4). Increasing the decoding radii means the decoder
% searches over a larger range of decoding set, so the performance will
% improve i.e. the error floors will drop with rl,ru.

% Recall rad_l, rad_u are non-negative integers so increment them by
% 1 at a time:
[floor_MD,floor_FA,floor_AUE] = ...
    RCU_floor_KaRandomUnknown_SRA([],rad_l,rad_u, tail_prob, n,L,E_Ka,p_Ka);
% Above means that for P,P1\to\infty, eps_MD/FA/AUE>=floor_MD/FA/AUE, 
% given the current rad_l,rad_u.

% If total error floor max(floor_MD,floor_FA)+floor_AUE > target_epsTotal,
% we need to increase rad_l,rad_u until the sign flips around:
while floor_MD+floor_AUE > target_epsTotal || floor_FA+floor_AUE > target_epsTotal 
    if floor_MD+floor_AUE > target_epsTotal
        rad_u = rad_u + 1; % because as stated in eq(38), floor_MD 
        % decreases by increasing rad_u
    end
    if floor_FA+floor_AUE > target_epsTotal
        rad_l = rad_l + 1; % because as stated in eq(40), for fixed rad_u, 
        % floor_FA decreases by increasing rad_l
    end
    [floor_MD,floor_FA,floor_AUE] = ...
        RCU_floor_KaRandomUnknown_SRA([],rad_l,rad_u,tail_prob,n,L,E_Ka,p_Ka);
end
% By now, the rad_l,rad_u ensure 
% max(floor_MD,floor_FA)+floor_AUE < target_epsTotal.

% To be conservative, increase search range even further if the target 
% is particularly small:
if target_epsTotal < 1e-1
    rad_l = rad_l + 2;
    rad_u = rad_u + 2;
end
fprintf(log_file, "Adjusted radii to rad_l=%e, rad_u=%e", rad_l, rad_u);
end