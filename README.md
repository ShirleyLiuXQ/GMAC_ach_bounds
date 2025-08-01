
# Finite-Length Random Coding Achievability Bounds for Many-User Random Access

This repository contains the MATLAB code that computes the finite-length achievability bounds in\
[1] X. Liu, P. Pascual Cobo, and R. Venkataramanan, "*Many-User Multiple Access with Random User Activity: Achievability Bounds and Efficient Schemes*," in submission, 2025.\
[2] X. Liu, P. Pascual Cobo, and R. Venkataramanan, "*Many-User Multiple Access with Random User Activity*," in Proc. IEEE Int. Symp. Inf. Theory, 2024.

This code is adapted from the source code by Khac-Hoang Ngo published at https://github.com/khachoang1412/UMA_random_user_activity.

---
We consider the real-valued Gaussian multiple access channels (GMAC), where the number of active users $K_a$ is random and unknown. This incurs three types of errors: misdetections (MD), false alarms (FA), and active user errors (AUE). This code computes these errors under a random coding scheme with maximum-likelihood decoding.

## Main functions are located in the folder `RCU_KaUnknown`:
- `RCU_KaRandomUnknown_SRA`: computes the finite-length achievability bounds on the three types of errors (Theorem 1 of [1]), allowing any distribution of $K_a$. The code computes the bounds for a given signal-to-noise ratio $E_b/N_0$ and decoding radius. "RCU" stands for "random coding union (bounds)", and "SRA" stands for "sourced random access".
- `RCU_KaBinomialUnknown_SRA`: computes the finite-length achievability bounds on the three types of errors (Theorem 1 of [1]), assuming $K_a$ is Binomially distributed. 
- `RCU_floor_KaRandomUnknown_SRA`: computes the error floors for the three types of errors (Corollary 1 of [1]), as $E_b/N_0$ tends to infinity.
- `EbN0_KaBinomialUnknown_SRA`: calculates the minimum $E_b/N_0$ needed for the achievability bounds to fall below a pre-specified target total error, defined as $\max(p_{\text{MD}},p_{\text{FA}})+p_{\text{AUE}}$.

Auxiliary functions: 
- `golden_search_P1_SRA`: searches for the optimal $P'$ in Theorem 1 of [1] to minimise the errors for a given $E_b/N_0$ and decoding radius.
- `binary_search_P_SRA`: searches for the smallest $P$ such that errors fall below a pre-specified threshold.

## Two example experiments:
`gmac_RCU`:  computes and plots the three types of errors for varying $E_b/N_0$.
`gmac_EbN0`: for a given active-user spectral efficiency, computes and plots the minimum $E_b/N_0$ needed for errors to fall below a pre-specified threshold.