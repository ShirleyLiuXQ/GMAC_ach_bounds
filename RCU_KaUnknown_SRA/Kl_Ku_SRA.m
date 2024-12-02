function [K_l, K_u] = Kl_Ku_SRA(L, E_Ka, p_Ka, tail_prob)
% Find K_l, K_u such that prob(Ka \notin [K_l, K_u]) is small.
% Ensure the same Kl, Ku are used for both computing the bounds and the
% asymptotic error floors.
% L    : total number of users
% E_Ka : expected number of active users Ka, E_Ka\in[0,L]
% p_Ka : PMF of Ka

% tail_prob = 1e-7; %1e-12
K_l = floor(E_Ka); 
K_u = ceil(E_Ka);
while p_Ka(K_l-1) > tail_prob
    K_l = K_l - 1;
end
K_l = max(K_l,0);
while p_Ka(K_u+1) > tail_prob
    K_u = K_u + 1;
end
assert(K_l>=0 && K_u<=L);
end