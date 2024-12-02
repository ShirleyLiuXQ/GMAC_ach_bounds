function data = RCU_KaBinomialUnknown_SRA(k, n, L, alpha_list, ...
    EbN0db_list, rad_l_list, rad_u_list, tail_prob, obj, ...
    P1_asFactorOfP, frac_TOL_golden)
% Compute the RCU bounds on pMD, pFA, pAUE for the L-user GMAC with the 
% number of active users unknown and following a Binomial distribution
% Bin(L,1-alpha) i.e. among the L users, each user is active with prob 
% 1-alpha.
%
% INPUTS
%   k   : number of bits per symbol
%   n   : framelength (number of real DoFs)
%   L   : total num users, among which Ka is active
%   alpha_list : set of probabilities that a user is silent, replacement of
%                E_Ka_list
%   EbN0db_list: energy per bit in dB for real signals
%   rad_l_list : vector storing set of values for the lower decoding radius
%   rad_u_list : vector storing the set of corresponding upper decoding radius
% Note: rad_l_list and rad_u_list are of the same length. The ith entries
% in the two vectors give the ith lower-upper-decoding-radius pair.
%   obj: can be 'max' or 'max_pTotal' which determines the objective 
%   function to minimise wrt P1.
%   P1_asfactorOfP: when empty, this function golden-searches for the P1 
%   that minimises obj. When non-empty, obj or golden-search is not used, 
%   this function directly evaluates the three types of errors using the 
%   optimal P1.
% 
%
% OUTPUTS
%   data : store the system parameters and the computed RCU bound i.e. the
%   minimum achievable eps_MD*, eps_FA*, eps_AUE* (and the corresponding
%   P1* that leads to eps_MD*, eps_FA*, eps_AUE*) for each given 
%   combination of E_Ka, EbN0db and the rad_l-rad_u pair.
% 
% This function requires EbN0 as input and returns eps_MD, eps_FA, eps_AUE
% for the given EbN0, while its counterpart EbN0_KaBinomialUnknown_SRA
% requires  target_eps_MD, target_eps_FA, target_eps_AUE and returns the
% minimal EbN0 needed to achieve these targets.

tStart = tic;
DEBUG = 0;

%% debugging mode
if DEBUG == 1
    k = 128;
    n = 19200;
    EbN0db_list = 4;
    L = 100;
    alpha_list = 0.5;
    % In debugging mode, set decoding radii to zero. This means the
    % minimisation in eq(5) is over W' whose size is exactly Ka'
    rad_l_list = 0;
    rad_u_list = 0;
end
assert(length(rad_l_list) == length(rad_u_list));
numRad = length(rad_l_list);

data.k      = k;
data.n      = n;
data.L      = L;
data.alpha  = alpha_list;
data.EbN0db = EbN0db_list; % increasing
data.rad_lower = rad_l_list;
data.rad_upper = rad_u_list;
data.pKa = 'Binomial';
data.P1_asFactorOfP = P1_asFactorOfP;
data.tail_prob = tail_prob;
if isempty(P1_asFactorOfP) % unknown the form of optimal P1
    data.obj = obj;
end

%% symbol power budget
% Eb/N0=nP/k when the signals are complex
% Eb/N0=nP/(2k) when the signals are real
P_list = 2*k.*10.^(EbN0db_list./10)./n; 
data.P = P_list;
numEbN0 = length(EbN0db_list);
if numEbN0 > 1
    assert(all(diff(EbN0db_list) >=0)); % monotonically increasing
end
numAlpha = length(alpha_list);

%% initialization
pMD = ones(numEbN0,numAlpha,numRad);
pFA = ones(numEbN0,numAlpha,numRad);
pAUE = ones(numEbN0,numAlpha,numRad);
P1 = zeros(numEbN0,numAlpha,numRad);
floor_pMD = ones(numAlpha,numRad);
floor_pFA = ones(numAlpha,numRad);
floor_pAUE = ones(numAlpha,numRad);

%% Compute the RCU bounds
for iAlpha = 1:numAlpha
    alpha = alpha_list(iAlpha);
    fprintf('alpha=%.2f [%d/%d]\n', alpha, iAlpha, numAlpha);
    p_Ka = @(K) binopdf(K,L,1-alpha); % 1-alpha active; alpha inactive
    E_Ka = L*(1-alpha);
    for iRad = 1:numRad
        rad_l = rad_l_list(iRad); 
        rad_u = rad_u_list(iRad); 
        fprintf('rad_l=%.2f,rad_u=%.2f [%d/%d]\n', rad_l, rad_u, iRad, numRad);
        [floor_pMDtmp,floor_pFAtmp,floor_pAUEtmp] = ...
            RCU_floor_KaRandomUnknown_SRA(P1_asFactorOfP,rad_l,rad_u, ...
            tail_prob, n,L,E_Ka,p_Ka);
        floor_pMD(iAlpha,iRad) = floor_pMDtmp;
        floor_pFA(iAlpha,iRad) = floor_pFAtmp;
        floor_pAUE(iAlpha,iRad) = floor_pAUEtmp;
        fprintf('\n[pMD_floor = %e]\n[pFA_floor = %e]\n[pAUE_floor = %e]\n', ...
            floor_pMDtmp, floor_pFAtmp, floor_pAUEtmp);

        data.floor_pMD = floor_pMD;
        data.floor_pFA = floor_pFA;
        data.floor_pAUE = floor_pAUE;
        for iEbN0 = 1:numEbN0
            fprintf('Eb/N0 [%d/%d]\n', iEbN0, numEbN0);
            P = P_list(iEbN0);
            if isempty(P1_asFactorOfP)
                fprintf('P1 is empty. Optimise over P1 to minimise the objective:');
                f_rcu = @(P,P1) RCU_KaRandomUnknown_SRA(P,P1, ...
                    rad_l,rad_u, tail_prob, k,n,L,E_Ka,p_Ka);
                % Optimise over P1 to minimise obj(p_MD,p_FA,pAUE):
                [pMDtmp,pFAtmp,pAUEtmp,P1tmp] = ...
                    golden_search_P1_SRA(f_rcu,P/2,P,P*frac_TOL_golden, obj,[],[],[]); 
                % start_val used to be 1e-9, changed to P/2 for speed
                % TOL used to be P/200. Changed to P/20 for speed.
                % Optimise over P1 to minimize a weighted sum of p_MD, p_FA, p_AUE:
                %     weight_MD = 1;
                %     weight_FA = 1;
                %     weight_AUE = 1;
                %     [p_MD_tmp,p_FA_tmp,p_AUE_tmp,P1_tmp] = ...
                %           golden_search_P1_SRA(f_rcu,0,P,P/100,'weighted', ...
                %           weight_MD,weight_FA,weight_AUE);
            else
                fprintf('P1 is specified. Skip optimisation over P1. Directly compute RCU bounds:')
                P1tmp = P1_asFactorOfP * P;
                [pMDtmp,pFAtmp,pAUEtmp] = RCU_KaRandomUnknown_SRA(...
                    P,P1tmp,rad_l,rad_u, tail_prob, k,n,L,E_Ka,p_Ka);
            end
            fprintf('\n[pMD = %e]\n[pFA = %e]\n[pAUE = %e]\n', pMDtmp, pFAtmp, pAUEtmp);
            pMD(iEbN0,iAlpha,iRad) = pMDtmp;
            pFA(iEbN0,iAlpha,iRad) = pFAtmp;
            pAUE(iEbN0,iAlpha,iRad) = pAUEtmp;
            P1(iEbN0,iAlpha,iRad) = P1tmp;
            
            if numEbN0 > 1
                % Early termination if pErrors have converged to their large-Eb/N0 asymptotics:
                [terminate, pMD(:,iAlpha,iRad), pFA(:,iAlpha,iRad), ...
                    pAUE(:,iAlpha,iRad)] = early_terminate(...
                    pMD(:,iAlpha,iRad), pFA(:,iAlpha,iRad), pAUE(:,iAlpha,iRad), iEbN0);
            end
            % Save results:
            sim_time = toc(tStart);
            data.pMD       = pMD;
            data.pFA       = pFA;
            data.pAUE      = pAUE;
            data.P1        = P1;
            data.sim_time  = sim_time;
            if numEbN0 > 1
                if terminate
                    break
                end
            end
        end
    end
end
pMD = squeeze(pMD);
pFA = squeeze(pFA);
pAUE = squeeze(pAUE);
P1 = squeeze(P1);
floor_pMD = squeeze(floor_pMD);
floor_pFA = squeeze(floor_pFA);
floor_pAUE = squeeze(floor_pAUE);

%% Save the results
sim_time = toc(tStart);
data.pMD       = pMD;
data.pFA       = pFA;
data.pAUE      = pAUE;
data.P1        = P1;
data.floor_pMD = floor_pMD;
data.floor_pFA = floor_pFA;
data.floor_pAUE = floor_pAUE;
data.sim_time  = sim_time;

% if DEBUG ~= 1
%     dt = datetime('now','TimeZone','local','Format','d-MMM-y_HH-mm-ss');
%     dtstr = char(dt); % string(dt) creates a string array; char(dt) a character array
%     filename = ['RCU_KaBinomialUnknown_SRA_obj_' obj '_alpha_' ...
%         num2str(min(alpha_list)) 'to' num2str(max(alpha_list)) ...
%         '_k_' num2str(k) '_n_' num2str(n) '_L_' num2str(L) ...
%         '_radL_' sprintf('%d', rad_l_list) ...
%         '_radU_' sprintf('%d', rad_u_list) ...
%         '_EbN0db_' num2str(min(EbN0db_list)) 'to' num2str(max(EbN0db_list)) '_' dtstr '.mat'];
%     save(filename, 'data', '-v7.3');
% else
%     keyboard
% end

end


function [terminate, pMD_arr, pFA_arr, pAUE_arr] = early_terminate(...
    pMD_arr, pFA_arr, pAUE_arr, iter_idx)
% Terminate when pErrors has converged to their large-Eb/N0 asymptotics.
% pMD_arr, pFA_arr, pAUE_arr need to be 1D arrays

len = length(pMD_arr);
assert((len == length(pFA_arr)) && (len == length(pAUE_arr)));
if iter_idx < 3
    terminate = false;
else
    assert((iter_idx >= 3) && (len >= iter_idx));
    atol = 0.05;
    approached_zero = max([pMD_arr(iter_idx), pFA_arr(iter_idx), pAUE_arr(iter_idx)] < atol);
    rtol = 0.05;
    ratio_MD = (pMD_arr(iter_idx-1) - pMD_arr(iter_idx)) / pMD_arr(iter_idx);
    ratio_FA = (pFA_arr(iter_idx-1) - pFA_arr(iter_idx)) / pFA_arr(iter_idx);
    ratio_AUE = (pAUE_arr(iter_idx-1) - pAUE_arr(iter_idx)) / pAUE_arr(iter_idx);
    has_converged = max([ratio_MD, ratio_FA, ratio_AUE]) < rtol;
    
    terminate = approached_zero && has_converged;
    % Fill in the remaining entries:
    pMD_arr(iter_idx+1: end) = pMD_arr(iter_idx);
    pFA_arr(iter_idx+1: end) = pFA_arr(iter_idx);
    pAUE_arr(iter_idx+1: end) = pAUE_arr(iter_idx);
end
end