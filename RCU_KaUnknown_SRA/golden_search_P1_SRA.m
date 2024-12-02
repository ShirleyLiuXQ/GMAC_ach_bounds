function [eps_MD, eps_FA, eps_AUE, P1] = golden_search_P1_SRA(f, ...
    start_val, end_val, TOL, obj_type, weight_MD, weight_FA, weight_AUE)
% Minimize a weighted sum of eps_MD, eps_FA, eps_AUE over P1 \in [0,P], 
% where eps_MD, eps_FA, eps_AUE are output from f(P,P1), which is the RCU
% bounds in our main theorem.
% 
% INPUTS
%   f: function handle that requires P,P1 as inputs and returns fMD, fFA, 
%      fAUE. f is a strictly unimodal function with an extremum inside the
%      interval [START_INT, END_INT]
%   START_INT, END_INT: search over P1\in[START_INT, END_INT] for the minimum 
%   TOL: stop the search when the search interval shrinks from 
%   [START_INT, END_INT] to an interval with size TOL
%   obj_type: type of objective to minimise. Can be 'max', 'weighted' or
%   'max_pTotal'.
%   - 'max': obj = max(pMD,pFA,pAUE); 
%   - 'weighted': obj = (weight_MD*pMD + weight_FA*pFA + weight_AUE*pAUE)/
%     (weight_MD + weight_FA + weight_AUE)
%   - 'max_pTotal': obj = max(pMD,pFA)+pAUE
%   weight_MD, weight_FA, weight_AUE: If type='weighted', specify the
%   weights on eps_MD, eps_FA, eps_AUE via weight_MD, weight_FA, weight_AUE
% 
% Adapted from golden_search which minimises f that takes in a single
% input and returns a single output. The only difference is that this
% function minimises a composite objective function which combines the 
% 3 quantities returned by f in a prescribed way. The 3 quantities all 
% depend on the same input variable P1, so the composite objective is also
% minimised wrt that input variable. 
% 
% OUTPUTS
%   eps_MD, eps_FA, eps_AUE: MD, FA, AUE probabilities corresponding to the
%   minimum composite objective function
%   P1: the optimal P1 that minimises the composite objective function

tic
fprintf(['Running golden_search_P1_SRA (i.e. search over P1 in[0,P]' ...
    'to minimise error objective)...\n']);
% Define the objective function to minimise:
if strcmp(obj_type,'max')
    g = @(fMD,fFA,fAUE) max([fMD,fFA,fAUE]);
elseif strcmp(obj_type,'weighted')
    assert(weight_MD>=0 && weight_FA>=0 && weight_AUE>=0 && ...
        (weight_MD+weight_FA+weight_AUE)>0);
    g = @(fMD,fFA,fAUE) (fMD*weight_MD + fFA*weight_FA + fAUE*weight_AUE)/...
        (weight_MD+weight_FA+weight_AUE); % my own way of defining weighted sum
elseif strcmp(obj_type,'max_pTotal')
    g = @(fMD,fFA,fAUE) max(fMD,fFA)+fAUE;
else
    error('Invalid obj_type\n');
end

P = end_val;                  % power budget
tau = double((sqrt(5)-1)/2);  % golden ratio, around 0.618
max_iter = 30; % ESIT used 30          % maximum number of iterations
i_iter = 0;                   % iteration index

x1 = start_val+(1-tau)*(end_val-start_val);            
x2 = start_val+tau*(end_val-start_val);
P1_list = [x1, x2];

[fMD_x1,fFA_x1,fAUE_x1] = f(P,x1); 
[fMD_x2,fFA_x2,fAUE_x2] = f(P,x2);
pMD_list = [fMD_x1, fMD_x2];
pFA_list = [fFA_x1, fFA_x2];
pAUE_list = [fAUE_x1, fAUE_x2];

while (abs(end_val-start_val)>TOL) && (i_iter<max_iter)
    fprintf('golden search num_iter=%d/max_iter=%d\n', i_iter, max_iter);
    if g(fMD_x1, fFA_x1, fAUE_x1) < g(fMD_x2, fFA_x2, fAUE_x2) % function
        % smaller at x1 than at x2
        end_val = x2; % shrink interval from above
        % Update x1, x2 to maintain the ratios between the four intervals:
        x2 = x1; 
        x1 = start_val+(1-tau)*(end_val-start_val); 
        % Update f(x1), f(x2):
        fMD_x2 = fMD_x1; 
        fFA_x2 = fFA_x1;
        fAUE_x2 = fAUE_x1;
        [fMD_x1,fFA_x1,fAUE_x1] = f(P,x1); 
        P1_list = [P1_list, x1];
        pMD_list = [pMD_list, fMD_x1];
        pFA_list = [pFA_list, fFA_x1];
        pAUE_list = [pAUE_list, fAUE_x1];
    else
        start_val = x1; % shrink interval from below
        % Update x1, x2:
        x1 = x2;
        x2 = start_val+tau*(end_val-start_val); 
        % Update f(x1), f(x2):
        fMD_x1 = fMD_x2;
        fFA_x1 = fFA_x2;   
        fAUE_x1 = fAUE_x2;
        [fMD_x2,fFA_x2,fAUE_x2] = f(P,x2);
        P1_list = [P1_list, x2];
        pMD_list = [pMD_list, fMD_x2];
        pFA_list = [pFA_list, fFA_x2];
        pAUE_list = [pAUE_list, fAUE_x2];
    end
    i_iter=i_iter+1;
end
fprintf('golden search converged: num_iter=%d/max_iter=%d\n', i_iter, max_iter);

% Return the P1, eps_MD, eps_FA, eps_AUE corresponding to the minimum 
% objective function, which is the weighted sum of eps_MD, eps_FA, eps_AUE:
if g(fMD_x1, fFA_x1, fAUE_x1) < g(fMD_x2, fFA_x2, fAUE_x2)
    P1 = x1;
    eps_MD = fMD_x1;
    eps_FA = fFA_x1;
    eps_AUE = fAUE_x1;
else
    P1 = x2;
    eps_MD = fMD_x2;
    eps_FA = fFA_x2;
    eps_AUE = fAUE_x2;
end
fprintf('[golden_search_P1_SRA completed in %.2f secs]\n', toc);
assert(g(eps_MD, eps_FA, eps_AUE) <= min(g(pMD_list, pFA_list, pAUE_list)), ...
    'Final espTotal should be smaller than previous evaluations');
end