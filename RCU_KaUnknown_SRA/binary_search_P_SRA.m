function [epsMD, epsFA, epsAUE, epsTotal, P,optP1, data] = ...
    binary_search_P_SRA(f,P_l,P_u,target_epsTotal,TOL_binary, ...
    P1_asFactorOfP, frac_TOL_golden)
% Search for P\in [P_l,P_u] such that the total error i.e. 
% max(eps_MD,eps_FA)+eps_AUEeps_MD \in [target_epsTotal-TOL, target_epsTotal]
% where eps_MD, eps_FA and eps_AUE are obtained from a f(P,P1).
% 
% The prerequisite is that epsTotal monotonically decreases with P.
% 
% P1_asFactorOfP: float or empty.
% For each P value, we can either minimise our objective function 
% max(eps_MD,eps_FA)+eps_AUE over P1 or choose P1 simply as a fraction of P
% (Note 0.96*P is found to be the optimal P1 in Ngo22's experiments, but 
% our optimal P1 is likely different and needs to be recalibrated).
% Binary search with fixed P1 runs significantly faster.
% 
% No maximum iteration limit for binary search.
% data: for debugging use.

assert(P_l>=0 & P_u>=0, "P_l,P_u must be non-negative");
tStart = tic;
fprintf('Running binary_search_P_SRA...\n');

% Store the following quantities for each iteration:
% The first two entries store the initialisation f(P_l), f(P_u):
P_list = [];
epsTotal_list = [];
optP1_list = [];

% find_eps_P1 is a function handle. It takes in P and returns the P1* which 
% minimises the specified objective function in terms of eps_MD, eps_FA 
% and epsAUE. It returns P1*, and the eps_MD*, eps_FA* and epsAUE*
% corresponding to P1*.
if isempty(P1_asFactorOfP)
    % We do not give flexibility in the form of the objective.
    find_eps_P1 = @(P) golden_search_P1_SRA(f,P/2,P, ...
        P*frac_TOL_golden,'max_pTotal',[],[],[]);
    % start_val used to be 1e-9, changed to P/2 for speed
    % TOL used to be P/200. Changed to P/20 for speed.
else
    find_eps_P1 = @(P) compute_err_fixedP1(f, P, P1_asFactorOfP);
end

%% Given fixed P1 or the optimal P1 we find, search for P that gives 
% epsilons below the targets:
[epsMD_l,epsFA_l,epsAUE_l, P1_l] = find_eps_P1(P_l);
epsTotal_l = max(epsMD_l, epsFA_l) + epsAUE_l;
% P_l is a lower bound on P (lower P leads to larger errors)
P_list=[P_list, P_l];
epsTotal_list=[epsTotal_list, epsTotal_l];
optP1_list=[optP1_list, P1_l];

[epsMD_u,epsFA_u,epsAUE_u, P1_u] = find_eps_P1(P_u); 
epsTotal_u = max(epsMD_u, epsFA_u) + epsAUE_u;
% x2 is a upper bound on P (higher P leads to smaller errors) 
P_list=[P_list, P_u];
epsTotal_list=[epsTotal_list, epsTotal_u];
optP1_list=[optP1_list, P1_u];
      
% The warnings in case 2 and 3 apply because find_eps_P1 is monotonic on
% P\in[P_l,P_u] in the total error probability max(eps_MD,eps_FA)+eps_AUE.
% Higher power P gives smaller total error, and vice versa.
i_iter = 0; 
if P_l == P_u || (epsTotal_u > target_epsTotal)
    if epsTotal_u > target_epsTotal
        warning('No power in range can achieve target error\n');
    end
    P = -1;
    epsMD = -1;
    epsFA = -1;
    epsAUE = -1;
    epsTotal = -1;
    optP1 = -1;
    data = save_data(i_iter, P_list, epsTotal_list, optP1_list);
elseif epsTotal_l < target_epsTotal
    warning('All power in range can achieve target error\n');
    P = P_l;
    epsMD = epsMD_l;
    epsFA = epsFA_l;
    epsAUE = epsAUE_l;
    epsTotal = epsTotal_l;
    optP1 = P1_l;
    data = save_data(i_iter, P_list, epsTotal_list, optP1_list);
else 
    % There exists a P*\in[P_l,P_u] such that epsTotal<target when P>P* and 
    % epsTotal may>target when P<P*. In this case, use binary search.
    % Use binary search until the target is found, without iteration limits.     
    while (epsTotal_l >= epsTotal_u)
        i_iter = i_iter+1;
        fprintf('binary search num_iter=%d\n', i_iter);
        x_mid = (P_l+P_u)/2;
        [epsMD_mid,epsFA_mid,epsAUE_mid,P1_mid] = find_eps_P1(x_mid); 
        epsTotal_mid = max(epsMD_mid, epsFA_mid) + epsAUE_mid;
        P_list=[P_list,x_mid];
        epsTotal_list=[epsTotal_list,epsTotal_mid];
        optP1_list=[optP1_list,P1_mid];
        % Recall that the goal is to find the minimal power that gives
        % total error below target:
        if (target_epsTotal >= epsTotal_mid && ...
                epsTotal_mid >= target_epsTotal-TOL_binary)
            P = x_mid;
            epsMD = epsMD_mid;
            epsFA = epsFA_mid;
            epsAUE = epsAUE_mid;
            epsTotal = epsTotal_mid;
            optP1 = P1_mid;
            fprintf('binary search converged num_iter=%d\n', i_iter);

            data = save_data(i_iter, P_list, epsTotal_list, optP1_list);
            fprintf('[binary_search_P_SRA completed in %.2f secs]\n', toc(tStart));
            return;
        elseif target_epsTotal > epsTotal_mid 
            P_u = x_mid; % target power is in the lower-half, so update UB
        else
            P_l = x_mid; % target power is in the upper-half, so update LB
        end
    end
    assert(false, 'binary search shouldnt reach here.');
end
end

function data = save_data(i_iter, P_list, epsTotal_list, optP1_list)
    data.num_iter_conv = i_iter; % convergence
     
    % Sort P_list in ascending order:
    [P_list, idx] = sort(P_list, 'ascend');
    data.P_list = P_list;
    data.epsTotal_list = epsTotal_list(idx);
    data.optP1_list = optP1_list(idx);
    epsTotal_decreases_w_P = all(diff(data.epsTotal_list) <=0);
    % issorted(data.epsTotal_list, 'descend')
    data.epsTotal_decreases_w_P = epsTotal_decreases_w_P;
    % expect epsTotal to be monotonic in P
    assert(data.epsTotal_decreases_w_P, "epsTotal should monotonically decrease with P\n");
end

function [eps_MD, eps_FA, eps_AUE, P1] = compute_err_fixedP1(f, P, fracP1)
% Wrapper function that returns the three types of errors and P1.
    P1 = P*fracP1; % Prescribed P1. No need for golden_search over P1.
    [eps_MD, eps_FA, eps_AUE] = f(P,P1); % The f we'll use is the f_rcu in 
    % EbN0_KaPoissonUnknown, which returns the three errors. 
end