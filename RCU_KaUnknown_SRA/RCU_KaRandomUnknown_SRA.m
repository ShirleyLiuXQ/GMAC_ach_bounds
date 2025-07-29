function [eps_MD,eps_FA,eps_AUE] = RCU_KaRandomUnknown_SRA(P,P1, ...
    rad_l,rad_u, tail_prob, k,n,L,E_Ka,p_Ka, log_file_name)
% Compute the RCU bounds on the misdetection, false alarm probabilities,
% and the active user error respectively, p_MD, p_FA, p_AUE, for the
% standard real-valued GMAC, where different users use different codebooks,
% with num active users random and unknown and with a fixed P1<P.
% Counterpart of RCU_KaRandomUnknown, using only error-exponent-type bounds
% ptt1, not the information-density bounds qtt1.
%
% [1]: Ngo22
%
% INPUTS
%   P   : symbol power budget
%   P1  : the actual power used in our random code (denoted by P' in Ngo22)
%   rad_l : lower decoding radius
%   rad_u : upper decoding radius
%   k     : number of bits per symbol
%   n     : framelength (number of real channel uses/ DoFs)
%   L     : total num users
%   E_Ka  : average number of active users
%   p_Ka  : PMF of the number of active users Ka
%
% OUTPUTS
%   eps_MD, eps_FA, eps_AUE : upper bounds on p_MD, p_FA, p_AUE given the
%   input parameters

log_file = fopen(log_file_name, 'a');
tStart = tic;
fprintf(log_file, ['Running RCU_KaRandomUnknown_SRA which calculates pErrors for ' ...
    'given P,P1,rad_l,rad_u...\n']);
% codebook size
M = 2^k;

% n real channel uses is equivalent to n/2 complex channel uses:
assert(mod(n,2) == 0, ['Num real channel uses n has to be even, but the ' ...
    'input n=%d is not even'], n);
n = n/2;

%% Computation of \tilde{p}
[K_l, K_u] = Kl_Ku_SRA(L, E_Ka, p_Ka, tail_prob);
ptilde = 1 - sum(p_Ka(K_l:K_u)) + E_Ka*gammainc(n*P/P1,n,'upper'); % our \tilde{p}

if ptilde >= 1
    fprintf(log_file, 'ptilde greater than 1, return\n')
    eps_MD = 1;
    eps_FA = 1;
    eps_AUE = 1;
    return
end

%% Initialize eps_MD, eps_FA, eps_AUE to be \tilde{p}
eps_MD = ptilde;
eps_FA = ptilde;
eps_AUE = ptilde;

%% The expectation w.r.t. Ka (outer sum over Ka):
fprintf(log_file, ['Starting parfor to compute summands parameterised by Ka ' ...
    'in range [Kl=%d, Ku=%d]\n'], K_l, K_u);
parfor Ka = K_l:K_u % parfor ignores the max{Kl,1} lower bound in eq(6)
    local_tStart = tic;
    fprintf('== Ka=%d ==', Ka);
    % Nevertheless the lower bound of 1 is enforced in sum_psi, by setting
    % the terms sumMD, sumAUE with Ka=0 on their denominators to zeros.
    %% Compute the term \xi(Ka,Ka') using [1, Th. 2]
    % The function \zeta for ML esimator of Ka (with K=Ka):
    KaEst_thres = @(Ka,KaEst,P1) n.*(log(1+Ka.*P1) - log(1+KaEst.*P1))./...
        ((1+Ka.*P1)./(1+KaEst.*P1)-1); % eq(29)

    % For energy-based estimator of Ka, use the below (eq(30))
    % KaEst_thres = @(Ka,KaEst,P1) n./(1+Ka.*P1).*(1+(Ka+KaEst).*P1/2);

    % Option 1: compute the exact term defined in [1]
    % tic
    % P_Ka_KaEst_list = zeros(Ka_u - Ka_l + 1,1);
    % for idx = 1:Ka_u - Ka_l + 1
    %     KaEst = Ka_l + idx - 1;
    %     P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,Ka_l:KaEst-1,P1),n,'lower');
    %     P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,KaEst+1:Ka_u,P1),n,'upper');
    %     P_Ka_KaEst_list(idx) = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b]);
    % end
    % toc

    % Option 2: slight relaxation that is faster
    % tic
    % 1) For all Ka'\neq Ka, instead of minimise over all K\in[Kl,Ku] as in
    % eq(28), (81), minimise subject to K=Ka. This gives a \xi(Ka,Ka')
    % thats suboptimal to the \xi in (28): (gammainc includes normalisation
    % factor Gamma(n))
    % suboptimal \xi(Ka,Ka') for Ka>Ka':
    P_Ka_KaEst_list_a = gammainc(KaEst_thres(Ka,K_l:Ka-1,P1),n,'lower');
    % suboptimal \xi(Ka,Ka') for Ka<Ka':
    P_Ka_KaEst_list_b = gammainc(KaEst_thres(Ka,Ka+1:K_u,P1),n,'upper');
    % P_Ka_KaEst_list_a, P_Ka_KaEst_list_b store the suboptimal estimates
    % of \xi(Ka,Ka') for the different Ka'\in[K_l,K_u] and Ka'\neq Ka.
    % 2) For Ka'=Ka, follow (28) exactly: see my notes p.16-17 for
    % justification of the expression below. The zero inserted to max
    % ensures max() is non-empty when P_Ka_KaEst_list_a P_Ka_KaEst_list_b
    % are both empty.
    P_Ka_Ka = 1 - max([P_Ka_KaEst_list_a P_Ka_KaEst_list_b 0]);
    % 3) P_Ka_KaEst_list stores \xi(Ka,Ka') (estimated or exact) for the
    % different Ka'\in[K_l,K_u]:
    P_Ka_KaEst_list = [P_Ka_KaEst_list_a'; P_Ka_Ka; P_Ka_KaEst_list_b'];
    % toc

    %% The expectation w.r.t. Ka' (inner sum over Ka'):
    for KaEst = K_l:K_u
        % Given Ka, Ka', \xi(Ka,Ka') and Kl, Ku are fixed.
        % \xi(Ka,Ka')
        xi_Ka_KaEst = P_Ka_KaEst_list(KaEst - K_l + 1);
        KaEst_l = max(KaEst - rad_l,K_l); % \underline{Ka'}
        KaEst_u = min(KaEst + rad_u,K_u); % \overline{Ka'}

        % Sets of possible values of t and t'
        t_vec = 0:min(KaEst_u,Ka); % our T
        num_t = length(t_vec);
        t1_l = max(Ka-KaEst_u,0) - max(Ka-KaEst_l,0); % 1st term in eq(23)
        t1_u = max(KaEst_u-Ka,0) - max(KaEst_l-Ka,0); % 1st term in eq(24)
        t1_mat = zeros(num_t, KaEst_u-KaEst_l+1);
        num_t1 = size(t1_mat,2);
        for i_t = 1:num_t
            % The set of values that t1 takes depends on t, so
            % we calculate and save the set of values of t1 for each t:
            t1_mat(i_t,:) = t1_l+t_vec(i_t):t1_u+t_vec(i_t);
            % Each row of t1_mat is \overline{T}_t in (23) without the zero
            % lower bound, and without 2nd UB in u_t.
            % Thus to use the exact def of \overline{T}_t in (23) means
            % we might need to trim entries off the ends of some rows in
            % t1_mat later on.
        end

        %% Computation of pt:

        % R_1 is R_1 in our original manuscript*t'. This makes the 1st term in E(t,t')
        % rho*rho1*R1 instead of rho*rho1*t'*R1.
        % R1_f = @(n,M,t1) t1/n * log(M); % used n corresponding to complex
        % Corrected R1_f during TIT revision: the new R1_f is defined to be
        % R1 in our original manuscript *t'.
        R1_f = @(n, M, L, Ka, KaEst_l, t, t1) calculate_R1(n, M, L, Ka, KaEst_l, t, t1); % used n corresponding to complex
        % eq(18):
        R2_f = @(n,Ka,K1,t) 1/n*(gammaln(min(Ka,K1)+1)-gammaln(t+1)-gammaln(min(Ka,K1)-t+1));

        % Trim off t1<0 entries:
        t1_mat2 = t1_mat;
        t1_mat2(t1_mat2<0) = 0;

        % First, consider t1 > 0, will later consider t1=0:
        % Initialize for the optimization over rho and rho1 \in[0,1]
        num_rho = 20; % ESIT used 50;
        rho_vec =  linspace(1e-9,1,num_rho);
        num_rho1 = 20; % ESIT used 50;
        rho1_vec = linspace(1e-9,1,num_rho1);
        % 4 params along 4 dims of the tensors: we have in order rho, rho1,
        % t, t1:
        rho = permute(repmat(rho_vec,length(rho1_vec),1,num_t,num_t1),[2,1,3,4]);
        rho1 = repmat(rho1_vec,length(rho_vec),1,num_t,num_t1);
        t = permute(repmat(t_vec,length(rho_vec),1,length(rho1_vec),num_t1),[1,3,2,4]);
        t1 = permute(repmat(t1_mat2,1,1,length(rho_vec),length(rho1_vec)),[3,4,1,2]);

        % Precompute some quantities: all using elementwise operations:
        P2 = 1+(max(Ka-KaEst_u,0)+max(KaEst_l-Ka,0)).*P1; % eq(16)
        P3 = P1.*(t1+rho.*t); % defined in remark iv) of Theorem 1
        P4 = rho.*rho1.*P3.*P2;
        P5 = rho1.*(rho-1).*t1.*P1;

        % Find optimal lambda, lambda*, that maximizes E0(t,t'). lambda* is
        % the root of a cubic equation as defined in remark iv) of Theorem 1:
        % c1_f...c4_f are as defined in remark iv):
        % (I assume c1_f...c4_f or Delta0_f...E0_f are correct - havent
        % checked through myself)
        c1_f = - P3.*P4.*P5 - (rho1+1).*t1.*P1.*P3.*P4;
        c2_f = P5.*(P3.^2 - P4) + rho1.*(t1.*P1.*P3.^2 - P3.*P4) ...
            - (2.*t1.*P1 + P3).*P4;
        c3_f = 2.*P3.*P5 + rho1.*(P3.^2 + t1.*P1.*P3) - 2.*P4;
        c4_f = P5 + rho1.*P3; % Checked

        % Below is the general cubic formula (see wikipedia Cubic equation
        % - General cubic formula:
        Delta0_f = c2_f.^2 - 3.*c1_f.*c3_f;
        Delta1_f = 2.*c2_f.^3 - 9.*c1_f.*c2_f.*c3_f + 27.*c1_f.^2.*c4_f;
        Q_f = ((Delta1_f + sqrt(Delta1_f.^2 - 4.*Delta0_f.^3))./2).^(1/3);
        lambda = real(-(c2_f + Q_f + Delta0_f./Q_f)./3./c1_f); % optimal lambda

        % Compute E0(rho,rho') as in eq(12): substituting lambda*
        E0_f = rho1.*(rho-1).*log(1+lambda.*P1.*t1) ...
            + (rho1-1).*log(1+lambda.*P3) ...
            + log(1+lambda.*P3 - lambda.^2.*P4);

        % Compute E(t,t') as in eq(11): maximise over rho, rho1 (the 1st
        % and 2nd dimensions).
        % Ett_f is a matrix of dim num_t*num_t1
        % Ett_f = reshape(max(max(-rho.*rho1.*R1_f(n,M,t1) ...
        %     - rho1.*R2_f(n,Ka,KaEst_u,t) + E0_f)), ...
        %     [num_t num_t1]); % Note the 1st term in E(t,t') is
        % rho*rho1*R1 instead of rho*rho1*t'*R1.
        % Corrected E(t,t') during TIT revision:
        Ett_f = reshape(max(max(-rho.*rho1.*R1_f(n,M,L,Ka,KaEst_l,t,t1) ...
            - rho1.*R2_f(n,Ka,KaEst_u,t) + E0_f)), ...
            [num_t num_t1]); 

        % Compute p_{t,t'} as in eq(10):
        % ptt is a num_t*num_t1 matrix. Along 1st dim t varies, along 2nd
        % dim t1 varies.
        ptt = min(exp(-n.*Ett_f),1); % cap each entry of the matrix by 1

        % Now, consider the case t1 = 0, where the optimal lambda is the
        % root of a quadratic function:
        rho = permute(repmat(rho_vec,length(rho1_vec),1,num_t),[2,1,3]);
        rho1 = repmat(rho1_vec,length(rho_vec),1,num_t);
        t = permute(repmat(t_vec,length(rho_vec),1,length(rho1_vec)),[1,3,2]);
        P3 = P1.*rho.*t; % def of P3 in Remark iv) after thm 1 with t1=0
        P4 = rho.*rho1.*P3.*P2;
        % The c2_f,c3_f,c4_f below are as defined in Remark iv) after thm 1
        % with t1=0:
        c2_f = -(rho1+1).*P3.*P4;
        c3_f = rho1.*P3.^2 - 2.*P4;
        c4_f = rho1.*P3;
        Delta = c3_f.^2 - 4.*c2_f.*c4_f; % b^2-4ac in quadratic formula
        lambda = -(c3_f+sqrt(Delta))./2./c2_f; % (-b+/-sqrt(D))/(2a)
        E0_f = (rho1-1).*log(1+lambda.*P3) + log(1+lambda.*P3 - lambda.^2.*P4);
        Ett_f = squeeze(max(max(- rho1.*R2_f(n,Ka,KaEst_u,t) + E0_f))); % a len-num_t vector
        ptt_t1zero = min(exp(-n.*Ett_f),1); % cap each entry of the matrix by 1

        % Recall ptt is num_t by num_t1, and t1_mat is also num_t by num_t1
        % with each row containing the values of t1 for a fixed t.
        % Truncate from above i.e. enforcing 2nd UB in (24) to carve out
        % our \overline{T}_t:
        ptt(t1_mat > KaEst_u-max(KaEst_l-Ka,0)) = 0;
        % by replacing the UB on each row in t1_mat with ut in (24)
        ptt(t1_mat < 0) = 0; % Truncate from below by imposing the zero cap in the LB of (23)
        % Now each row of ptt stores the nonzero ptt's for a given t,
        % across t1\in\overline{T}_t, and stores zeros everywhere else.
        assert(sum(isnan(ptt), 'all') == 0);
        assert(sum(ptt>1, 'all') == 0);

        % Combine the cases t1 > 0 and t1 = 0:
        for i_t = 1:num_t
            for i_t1 = 1:num_t1
                if t1_mat(i_t,i_t1) == 0 % if t1=0 (cannot use t1_mat2
                    % because some zeros in t1_mat2 actually represent NAs)
                    % Replace the solution of the cubic equation with that
                    % of the quadratic, which should be more accurate:
                    ptt(i_t,i_t1) = ptt_t1zero(i_t);
                end
            end
        end

        % Compute pt: by now we've made sure each row of ptt contains the
        % p_{t,t'} for a fixed t and t'\in\overline{T}_t.
        % Summing along 2nd dim (i.e. along which t1 varies) gives pt for
        % different t values. Then cap each pt by 1.
        %         pt = min(sum(ptt,2),1); % a len-num_t vector

        % Take the min of pt, qt, and xi(Ka,Ka') for given Ka,Ka':
        % (pt here is already the min of pt and qt)
        %         pt = min(pt,xi_Ka_KaEst); % len-num_t vector
        ptt = min(ptt,xi_Ka_KaEst); % num_t*num_t1 matrix
        assert(sum(isnan(ptt), 'all') == 0);
        assert(sum(ptt>1, 'all') == 0);

        %% pmf over psi, new to the standard GMAC setting:
        t_mat = repmat(t_vec,num_t1,1).'; % num_t*num_t1 matrix
        % sumMD, sumFA, sumAUE are num_t*num_t1 matrices (summands over psi):
        [sumMD, sumFA, sumAUE] = sum_psi(L,M,Ka,KaEst_l,KaEst_u,t_mat,t1_mat,t1_mat2);

        %% Compute epsilon_MD and epsilon_FA for a given P'< P:
        if Ka > 0 % if Ka<=0, the current (Ka,Ka')-pair doesnt contribute to epsMD and epsAUE.
            % Add a summand (which itself is a sum over t,t1,psi):
            eps_MD = eps_MD + feval(p_Ka,Ka)*sum(sum(ptt.*sumMD, 2)); % sum over t1 then t
            eps_AUE = eps_AUE + feval(p_Ka,Ka)*sum(sum(ptt.*sumAUE), 2);
        end
        eps_FA = eps_FA + feval(p_Ka,Ka)*sum(sum(ptt.*sumFA, 2));
    end
    %     if eps_MD >= 1 && eps_FA >= 1
    %         break;
    %     end
fprintf('[%f secs]\n', toc(local_tStart));
end
eps_MD = min(eps_MD,1);
eps_FA = min(eps_FA,1);
eps_AUE = min(eps_AUE,1);
fprintf(log_file, 'epsMD=%e, epsFA=%e, epsAUE=%e\n', eps_MD, eps_FA, eps_AUE);
fprintf(log_file, '[RCU_KaRandomUnknown_SRA finished in %f secs]\n', toc(tStart));
end


function R1 = calculate_R1(n, M, L, Ka, KaEst_l, t, t1)
% Calculate tmin, r, and R1(t, t1) where R1(t, t1) is the corrected term in
% our TIT revision.
% Input:
% - L, Ka, KaEst_l, M, n: scalars
% - n: number of complex channel uses
% - t, t1: num_rho*num_rho1*num_t*num_t1 4D tensor. they are invariant
% across the rho and rho1 axes. There may exist some non-valid t1 values in
% the t1 input tensor, but thats okay, lines 225 and 227 will discard 
% corresponding entries.
% Output:
% - R1: num_rho*num_rho1*num_t*num_t1 4D tensor (tmin, r are also num_rho*num_rho1*num_t*num_t1 4D tensors)

assert(length(n)==1 && length(M)==1 && length(L)==1 && length(Ka)==1 && ...
    length(KaEst_l)==1, "n, M, L, Ka, KaEst_L are scalars");
assert(isequal(size(t), size(t1)), "t and t1 have equal size");
% Calculate pairwise minimum
tmin = ((t + t1) - abs(t - t1)) / 2;
% Calculate r
r = L - Ka + tmin - max(KaEst_l - Ka, 0) - max(t1 - t, 0);
% Calculate R1
R1 = (1/n) * (tmin * log(M) + gammaln(r+1) - gammaln(tmin+1) - gammaln(r-tmin+1));
assert(isequal(size(t), size(R1)), "t and R1 have equal size");
end

function [sumMD, sumFA, sumAUE] = sum_psi(L,M,Ka,KaEst_l,KaEst_u,t,t1,t1_nonneg)
% Input:
%   - L,M,Ka,KaEst_l,KaEst_u: scalars
%   - t,t1,t1_nonneg: num_t*num_t1 2D matrices (t1 may contain redundant
%     entries outside of \overline{T}_t; while t1_nonneg is equivalent to
%     t1 with zero cap applied from below)
% Output:
%   sumMD, sumFA, sumAUE: num_t*num_t1 matrices storing the summands in
%   epsMD, epsFA, epsAUE that have summed over t,t1,psi, for a given
%   (Ka,Ka')-pair.
%
% L,M,Ka,KaEst_l are all scalars:
assert(length(L)==1 && length(M)==1 && length(Ka)==1 && ...
    length(KaEst_l)==1 && length(KaEst_u)==1);
% Ensure t,t1 are matrices of the same dims:
assert(ismatrix(t) && isequal(size(t), size(t1_nonneg)) ...
    && isequal(size(t), size(t1)));
assert(sum(t<0, 'all') == 0 && sum(t1_nonneg<0, 'all') == 0);
[num_t, num_t1] = size(t);
tmin = ((t+t1_nonneg) - abs(t-t1_nonneg))/2; % entrywise minimum, psi\in[0,tmin]
max_tmin = max(max(tmin));
max_num_psi = max_tmin+1;
t = repmat(t,1,1,max_num_psi);
t1_nonneg = repmat(t1_nonneg,1,1,max_num_psi);
t1 = repmat(t1,1,1,max_num_psi);
tmin = repmat(tmin,1,1,max_num_psi);
psi = permute(repmat(0:max_tmin, num_t,1,num_t1), [1,3,2]);
if max_num_psi > 1
    assert(isequal(size(psi), [num_t,num_t1,max_num_psi]));
else
    assert(isequal(size(psi), [num_t,num_t1]));
end

% Truncate excessive psi entries along the 3rd dim to ensure the downstream
% nchoosek's dont receive negative inputs:
psi1 = psi;
psi(psi>tmin) = 0; % entrywise comparison
assert(sum(psi(tmin==0) ~=0) == 0);

r = L-Ka+tmin-max(KaEst_l-Ka,0)-max(t1_nonneg-t,0);
assert(sum(r<0, 'all') == 0);

% nu_norm_tmp = nu_f(M,Ka,KaEst_l,KaEst_u, r, t1, tmin, psi, psi1);
nu_norm = nu_f_stable(M,Ka,KaEst_l,KaEst_u, r, t1, tmin, psi, psi1);

% num_nans_tmp = sum(isnan(nu_norm_tmp), 'all');
% num_nans = sum(isnan(nu_norm), 'all');
% if num_nans_tmp == 0 && num_nans == 0
%     assert(max(abs(nu_norm_tmp-nu_norm), [], 'all') <1e-8);
% else
%     assert(num_nans_tmp > 0);
%     assert(num_nans == 0);
% end
%% Summands over psi in eps_MD, eps_FA, eps_AUE
if Ka>0 % if Ka<=0, the current (Ka,Ka')-pair doesnt contribute to epsMD or epsAUE.
    sumMD = sum(nu_norm.*(max(Ka-KaEst_u,0) + max(t-t1_nonneg,0) + psi)/Ka, 3);
    sumAUE = sum(nu_norm.*(tmin-psi)/Ka, 3);
else
    sumMD = zeros(num_t, num_t1);
    sumAUE = zeros(num_t, num_t1);
end
Mrx = Ka-t+t1_nonneg + max(KaEst_l-Ka,0) - max(Ka-KaEst_u,0); % \hat{Ka}
assert(sum(Mrx < KaEst_l | Mrx > KaEst_u, 'all') == 0);
% \hat{Ka}\in[KaEst_l:KaEst_u]\in[K_l:K_u]\in[0:L]
Mrx(Mrx == 0) = inf; % when \hat{Ka}=0, num FAs=0, so current (Ka,Ka')-pair
% makes no contribution to epsFA.
sumFA = sum(nu_norm.*(max(KaEst_l-Ka,0) + max(t1_nonneg-t,0) + psi)./Mrx, 3);
assert(sum(sumFA(Mrx(:,:,1) == 0) ~=0) == 0);
assert(isequal(size(sumMD), [num_t, num_t1]) && ...
    isequal(size(sumFA), [num_t, num_t1]) && ...
    isequal(size(sumAUE), [num_t, num_t1]));
assert(sum(sumMD<0, 'all') == 0 && ...
    sum(sumFA<0, 'all') == 0 && ...
    sum(sumAUE<0, 'all') == 0);
assert(sum(isnan(sumMD), 'all') == 0 && sum(isnan(sumFA), 'all') == 0 && ...
    sum(isnan(sumAUE), 'all') == 0);
end

function nu_norm = nu_f_stable(M,Ka,KaEst_l,KaEst_u, r, t1, tmin, psi, psi1)
% t1: may contain negative entries, while t1_nonneg doesnt
% psi1 may contain entries greater that the respective tmin entries, while
% psi doesnt
% r: number remaining codebooks, not the decoding radius
% nu_norm: nu(tmin,psi) normalised by sum_psi nu(tmin,psi)
assert(isequal(size(tmin), size(psi)));
if ndims(tmin) == 3
    [num_t,num_t1,max_num_psi] = size(tmin);
else
    [num_t,num_t1] = size(tmin);
    max_num_psi = 1;
end
    function ln_nCk = ln_nCk_f(n,k)
        % Use gammaln(A) to avoid underflow/ overflow in log(gamma(A)):
        assert(isequal(size(n), size(k)));
        ln_nCk_nGreater_f = @(n,k) gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1);
        ln_nCk = zeros(size(n));
        nIsGreater = n>=k;
        ln_nCk(nIsGreater) = ln_nCk_nGreater_f(n(nIsGreater), k(nIsGreater));
        ln_nCk(~nIsGreater) = -Inf;
    end

ln_nu = ln_nCk_f(r-tmin,psi) + psi*log(M) + ln_nCk_f(tmin,tmin-psi) + (tmin-psi)*log(M-1);
if max_num_psi > 1
    assert(isequal(size(ln_nu), [num_t,num_t1,max_num_psi]));
else
    assert(isequal(size(ln_nu), [num_t,num_t1]));
end
assert(sum(isnan(ln_nu), 'all') == 0);

% The entries in nu corresponding to tmin=0 are 1:
assert(sum(ln_nu(tmin==0) ~= 0) == 0);
% The entries in nu corresponding to r-tmin<psi are zeros:
assert(sum(ln_nu(r-tmin<psi) ~= -Inf) == 0);

% Set redundant nu(tmin,psi) entries to zero along the 3rd dim:
% Entrywise comparison:
ln_nu(psi1>tmin) = -Inf;
ln_nu(t1<0) = -Inf; % t is guaranteed to be >=0 but t1 isnt
ln_nu(t1 > KaEst_u-max(KaEst_l-Ka,0)) = -Inf;
assert(sum(isnan(ln_nu), 'all') == 0);

% Subtract common exponent from exp(a)/sum_iexp(i):
max_ln_nu = max(ln_nu, [], 3); % max along 3rd dim i.e. across different psi
assert(isequal(size(max_ln_nu), [num_t,num_t1]));
max_ln_nu(max_ln_nu == -Inf) = 0; % avoid operations on inf e.g. inf-inf,
% which gives nan
% max_ln_nu(max_ln_nu == inf) = 0;
ln_nu_sub_max = ln_nu - max_ln_nu;
assert(sum(ln_nu_sub_max>0, 'all') == 0);
nu_div_max = exp(ln_nu_sub_max);
assert(sum(isnan(nu_div_max), 'all') == 0);

% Calculate nu(tmin), the normalisation factor in the denominator:
nu_den_div_max = sum(nu_div_max, 3);
assert(isequal(size(nu_den_div_max), [num_t,num_t1]));

assert(sum(nu_den_div_max(t1(:,:,1) >= 0) == 0) == 0); % denominator
% entries for valid t1 shouldnt be zeros
nu_den_div_max(t1(:,:,1) < 0) = Inf; % denominator entries for invalid t1
% should be Inf so that the respective summand=0
nu_den_div_max(nu_den_div_max==0) = Inf; % respective summand=0 as they
% shouldnt exist

nu_norm = nu_div_max./nu_den_div_max;
if max_num_psi>1
    assert(isequal(size(nu_norm), [num_t,num_t1,max_num_psi]));
else
    assert(isequal(size(nu_norm), [num_t,num_t1]));
end
assert(sum(isnan(nu_norm), 'all') == 0);
assert(sum((nu_norm<0|nu_norm>1), 'all') == 0);
boolTmp = (tmin(:,:,1) == 0) + (t1(:,:,1) >= 0) == 2; % tmin=0 and t1 valid
assert(sum(nu_norm(boolTmp) ~= 1) == 0);
assert(sum(nu_norm(t1 < 0) ~= 0) == 0);
assert(sum(nu_norm(t1 > KaEst_u-max(KaEst_l-Ka,0))~=0) == 0);
assert(sum(nu_norm(psi1>tmin) ~= 0) == 0);
assert(sum(nu_norm(r-tmin<psi) ~= 0) == 0);
end

function nu_norm = nu_f(M,Ka,KaEst_l,KaEst_u, r, t1, tmin, psi, psi1)
% nchoosek_f allows tensor inputs, unlike Matlab's builtin nchoosek:
% (Ensure n,k are tensors of same dims. Note nchoosek_f(0,0)=1,
% nchoosek_f(>0,0)=1, nchoosek_f(<0,any)=nan, nchoosek_f(>0,<0)=0,
% nchoosek_f(a,b>a)=0)
assert(isequal(size(tmin), size(psi)));
if ndims(tmin) == 3
    [num_t,num_t1,max_num_psi] = size(tmin);
else
    [num_t,num_t1] = size(tmin);
    max_num_psi = 1;
end

nCk_f = @ (n,k) gamma(n+1)./gamma(k+1)./gamma(n-k+1);

nu = nCk_f(r-tmin, psi).* M.^psi.*nCk_f(tmin,tmin-psi).*(M-1).^(tmin-psi);
if max_num_psi > 1
    assert(isequal(size(nu), [num_t,num_t1,max_num_psi]));
else
    assert(isequal(size(nu), [num_t,num_t1]));
end
try
    assert(sum(isnan(nu), 'all') == 0);
catch
    %     fprintf('nu calculated by nu_f contains nans; return nu_norm=all nan');
    nu_norm = nan;
    return
end
% The entries in nu corresponding to tmin=0 are 1:
assert(sum(nu(tmin == 0) ~= 1) == 0);
assert(sum(nu(r-tmin<psi) ~= 0) == 0);
% nu is a 3D tensor with num_psi pages, where each page is nu(tmin,psi) for
% the num_t*num_t1 different tmin values calculated using t, t1_nonneg, and
% a fixed psi.
% Set redundant nu(tmin,psi) entries to zero along the 3rd dim:
% Entrywise comparison:
nu(psi1>tmin) = 0;
nu(t1<0) = 0; % t is guaranteed to be >=0 but t1 isnt
nu(t1 > KaEst_u-max(KaEst_l-Ka,0)) = 0;
% We did not truncate nu(t1 > KaEst_u-max(KaEst_l-Ka,0)) = 0; according to
% the 2nd UB in ut eqn (24), because this has been taken care of by
% truncating ptt. Even if nu contains nonzero entries where
% t1 > KaEst_u-max(KaEst_l-Ka,0), their respective entries in ptt will be
% zeros, so they wont end up contributing to the sums.
assert(sum(isnan(nu), 'all') == 0);

% Calculate nu(tmin), the normalisation factor in the denominator:
nu_den = sum(nu, 3); % num_t*num_t1 matrix
boolTmp = (tmin(:,:,1) == 0) + (t1(:,:,1) >= 0) == 2; % t1_nonneg is t1
% with negative entries replaced by zeros
assert(sum(nu_den(boolTmp) ~= 1) == 0);
assert(sum(nu_den(t1(:,:,1) < 0) ~= 0) == 0);
assert(sum(nu_den(t1(:,:,1) >= 0) == 0) == 0);
nu_den(nu_den==0) = inf; % make the corresponding summand (which shouldnt
% exist)=0

% For large M: use Vandermondes convolution identity to approximate
% nu(tmin): (constants are numerically determined)
% if M > 2^10
%     nu_den1 = nchoosek_f(r(:,:,1),t(:,:,1)).*M.^t(:,:,1);
%     assert(max(max(abs(nu_den1-nu_den)./ nu_den1)) < 1e-3, ...
%         ['nu_den approximated by Vandermondes convolution identity ' ...
%         'is >1e-3 different from the exact nu_den']);
% end

nu_norm = nu./nu_den; % normalised nu
assert(sum(isnan(nu_norm), 'all') == 0);
assert(sum((nu_norm<0|nu_norm>1), 'all') == 0);
assert(sum(nu_norm(t1 > KaEst_u-max(KaEst_l-Ka,0))~=0) == 0);
assert(sum(nu_norm(t1 < 0) ~= 0) == 0);
assert(sum(nu_norm(psi1>tmin) ~= 0) == 0);
assert(sum(nu_norm(r-tmin<psi) ~= 0) == 0);
end
