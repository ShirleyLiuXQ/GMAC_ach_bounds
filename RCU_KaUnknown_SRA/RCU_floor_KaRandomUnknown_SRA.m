function [floor_MD,floor_FA,floor_AUE] = ...
    RCU_floor_KaRandomUnknown_SRA(P1_asFactorOfP,rad_l,rad_u, tail_prob, ...
    n,L,E_Ka,p_Ka)
% Compute the asymptotic (P, P1\to \infty) error floors for the RCU bounds 
% of the MD, FA and AUE probabilities for a system with random and unknown 
% number of active users. These floors are obtained by assuming i) P,P'\to
% \infty and ii) setting t,t'=0. 
% The expressions for floor_MD,floor_FA are given in Corollary 1 of Ngo22 
% and apply to our setting with adjusted \bar{p}. 
% Note eps_AUE=pbar when t,t'=0 so  floor_AUE is simply \bar{p}. 
%
% INPUTS
%   rad_l : lower decoding radius, in the definition of KaEst_l below eq(5)
%   rad_u : upper decoding radius, in the definition of KaEst_u below eq(5)
%   n     : framelength (number of real DoFs/ channel uses)
%   E_Ka  : average number of active users
%   p_Ka  : PMF of the number of active users Ka, a function handle that
%           admits one scalar/ vector input and returns the pmf evaluated
%           at its input.
%
% OUTPUTS
%   floor_MD, floor_FA, floor_AUE : the asymptotic error floors, that is,
%   for P, P1\to \infty, eps_MD>=floor_MD, eps_FA>=floor_FA and
%   eps_AUE>=floor_AUE.
% 
% Compared to RCU_KaRandomUnknown_SRA, this function requires two fewer
% inputs (i.e. P and P1) because the error floors are asymptotic, when P 
% and P1\to\infty.

tStart = tic;
fprintf('Running RCU_floor_KaRandomUnknown_SRA ...');

assert(mod(n,2) == 0, ['Num real channel uses n has to be even, but the ' ...
    'input n=%d is not even'], n);
n = n/2; % n real channel uses = n/2 complex channel uses

%% Calculate \bar{p}:
[K_l, K_u] = Kl_Ku_SRA(L, E_Ka, p_Ka, tail_prob);
pbar = 1 - sum(p_Ka(K_l:K_u)); % our \bar{p}
if ~isempty(P1_asFactorOfP)
    pbar = pbar + E_Ka*gammainc(n*1/P1_asFactorOfP,n,'upper');
%     pbar = pbar + E_Ka*gammainc(n*P/P1,n,'upper');
end
%% Initialize the floors to be \bar{p}:
floor_MD = pbar;
floor_FA = pbar;
floor_AUE = pbar;

if pbar < 1
    %% The outer sum over Ka in eq(31),(32):
    parfor Ka = K_l:K_u
        % Ka is the real num active users; KaEst (Ka') is the initial estimate
        % of Ka from y based on eq(4):
        % Compute \xi(Ka,Ka'):
        KaEst_thres = @(Ka,KaEst) n.*log(Ka./KaEst)./(Ka./KaEst-1); % this is
        % \zeta(K=Ka,Ka,Ka') for ML estimate of Ka in Thm 3
        % The justification for the following code is the same as its
        % counterpart in RCU_KaRandomUnknown:
        % suboptimal \xi(Ka,Ka') for Ka>Ka':
        P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,K_l:Ka-1),n,'lower'); 
        % suboptimal \xi(Ka,Ka') for Ka<Ka':
        P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,Ka+1:K_u),n,'upper');  
        % P_Ka_KaEst_list_a, P_Ka_KaEst_list_b are both >0, the zero inserted
        % in max ensures max() is non-empty
        P_Ka_Ka = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b 0]);
        % P_Ka_KaEst_list stores \xi(Ka,Ka') (estimated or exact) for the
        % different Ka'\in[K_l,K_u]:
        P_Ka_KaEst_list = [P_Ka_KaEst_list_a'; P_Ka_Ka; P_Ka_KaEst_list_b'];
    
        %% The inner sum over Ka' in eq(31),(32):
        for KaEst = K_l:K_u
            % P_Ka_KaEst is \xi(Ka, Ka'):
            if KaEst == 0
                P_Ka_KaEst = 0;
            else
                P_Ka_KaEst = P_Ka_KaEst_list(KaEst - K_l + 1);
            end
            % KaEst_l and KaEst_u are defined just below eq(5):
            KaEst_l = max(KaEst - rad_l, K_l); % \underline{K_a'}
            KaEst_u = min(KaEst + rad_u, K_u); % \overline{K_a'}
    
            if Ka > 0 % Ka=0 contributes a zero summand to floor_MD
                floor_MD = floor_MD + feval(p_Ka,Ka)*max(Ka-KaEst_u,0)*P_Ka_KaEst/Ka;
            end
            % The denominator in the inner sum over Ka' in (32):
            Mrx = (Ka + max(KaEst_l-Ka,0) - max(Ka-KaEst_u,0)); % number of decoded codewords
            if Mrx > 0 % Mrx=\hat{Ka}=0 contributes a zero summand to floor_FA
                floor_FA = floor_FA + feval(p_Ka,Ka)*max(KaEst_l-Ka,0)*P_Ka_KaEst/Mrx;
            end
        end
    end
end
floor_MD = min(floor_MD,1);
floor_FA = min(floor_FA,1);
floor_AUE = min(floor_AUE,1);
fprintf('\n[pMD_floor = %e]\n[pFA_floor = %e]\n[pAUE_floor = %e]\n', ...
            floor_MD, floor_FA, floor_AUE);

fprintf('[RCU_floor_KaRandomUnknown_SRA finished in %f secs]\n', toc(tStart));
end