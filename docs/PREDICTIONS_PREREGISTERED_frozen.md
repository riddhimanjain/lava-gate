# PRE-REGISTERED PREDICTIONS — gate-induced bias in LAVA local rg

**Generated:** 2026-08-14 09:03 UTC  
**Source:** closed form in `01_analytic.py`, no simulation involved.

> Frozen before `03_tierA_sim.py` was written or run. If the simulation disagrees with these numbers, the disagreement is the result and gets reported as such. This file is not to be edited afterwards.


## Parameters

- Reference panel N = 100,000 (from LD panel `.info` NOBS ~ 99,339)
- Univariate gate: chi-square, alpha = 0.05 (both traits must pass)
- theta_i = Omega_ii/Sigma_ii, the local signal-to-noise ratio
- param.lim = 1.25, cap.estimates = TRUE (applied in simulation, not here)


## Falsifiable claims

**P1 — Direction.** E[rho_hat | G] is attenuated toward zero: |E[rho_hat|G]| < |rho| for every scenario with r_e = 0.

**P2 — Mechanism is the denominator.** Attenuation equals 1/sqrt(a_1*a_2) to leading order, where a_i = E[omega_hat_ii|G]/Omega_ii is the conditional inflation of local h2. The numerator plays only a second-order role.

**P3 — Near-multiplicativity.** Because a_i depends on Omega_ii but not on Omega_12, the attenuation factor is (to leading order) independent of rho. So E[rho_hat|G] should be close to LINEAR through the origin in rho, with residuals from a straight-line fit below 0.01 in absolute value when r_e = 0.

**P4 — Sample overlap creates an intercept.** With r_e > 0, E[rho_hat|G] != 0 at rho = 0. The intercept grows with r_e.

**P5 — K dependence.** Attenuation worsens (slope falls) as K grows at fixed theta, because the gate threshold tightens relative to the signal.


## Predicted slope and intercept of E[rho_hat|G] against rho

| Scenario | K | theta1 | theta2 | r_e | a1 | a2 | 1/sqrt(a1 a2) | pred slope | pred intercept | max resid |
|---|---|---|---|---|---|---|---|---|---|---|
| S1_asdlike_x_sczlike | 300 | 0.065 | 0.22 | 0.0 | 2.961 | 1.161 | 0.5392 | 0.5897 | 0.0000 | 0.00000 |
| S2_matched_moderate | 300 | 0.1 | 0.1 | 0.0 | 2.023 | 2.023 | 0.4944 | 0.5260 | 0.0000 | 0.00000 |
| S3_overlap_re05 | 300 | 0.065 | 0.22 | 0.05 | 2.961 | 1.161 | 0.5392 | 0.5897 | 0.0001 | 0.00000 |
| S4_overlap_re10 | 300 | 0.065 | 0.22 | 0.1 | 2.961 | 1.161 | 0.5392 | 0.5897 | 0.0002 | 0.00000 |
| S5_smallK | 100 | 0.065 | 0.22 | 0.0 | 5.072 | 1.716 | 0.3390 | 0.5115 | 0.0000 | 0.00000 |
| S6_largeK | 500 | 0.065 | 0.22 | 0.0 | 2.360 | 1.051 | 0.6351 | 0.6605 | 0.0000 | 0.00000 |

## Predicted gate pass rates

| Scenario | P(trait1 passes) | P(trait2 passes) | P(both, indep approx) |
|---|---|---|---|
| S1_asdlike_x_sczlike | 0.1875 | 0.7862 | 0.1474 |
| S2_matched_moderate | 0.3141 | 0.3141 | 0.0987 |
| S3_overlap_re05 | 0.1875 | 0.7862 | 0.1474 |
| S4_overlap_re10 | 0.1875 | 0.7862 | 0.1474 |
| S5_smallK | 0.1189 | 0.4248 | 0.0505 |
| S6_largeK | 0.2397 | 0.9263 | 0.2220 |

## Full predicted curve

| Scenario | rho | predicted E[rho_hat \| G] |
|---|---|---|
| S1_asdlike_x_sczlike | 0.0 | 0.0000 |
| S1_asdlike_x_sczlike | 0.1 | 0.0590 |
| S1_asdlike_x_sczlike | 0.2 | 0.1179 |
| S1_asdlike_x_sczlike | 0.3 | 0.1769 |
| S1_asdlike_x_sczlike | 0.4 | 0.2359 |
| S1_asdlike_x_sczlike | 0.5 | 0.2948 |
| S1_asdlike_x_sczlike | 0.6 | 0.3538 |
| S1_asdlike_x_sczlike | 0.7 | 0.4128 |
| S1_asdlike_x_sczlike | 0.8 | 0.4717 |
| S1_asdlike_x_sczlike | 0.9 | 0.5307 |
| S2_matched_moderate | 0.0 | 0.0000 |
| S2_matched_moderate | 0.1 | 0.0526 |
| S2_matched_moderate | 0.2 | 0.1052 |
| S2_matched_moderate | 0.3 | 0.1578 |
| S2_matched_moderate | 0.4 | 0.2104 |
| S2_matched_moderate | 0.5 | 0.2630 |
| S2_matched_moderate | 0.6 | 0.3156 |
| S2_matched_moderate | 0.7 | 0.3682 |
| S2_matched_moderate | 0.8 | 0.4208 |
| S2_matched_moderate | 0.9 | 0.4734 |
| S3_overlap_re05 | 0.0 | 0.0001 |
| S3_overlap_re05 | 0.1 | 0.0591 |
| S3_overlap_re05 | 0.2 | 0.1181 |
| S3_overlap_re05 | 0.3 | 0.1770 |
| S3_overlap_re05 | 0.4 | 0.2360 |
| S3_overlap_re05 | 0.5 | 0.2950 |
| S3_overlap_re05 | 0.6 | 0.3539 |
| S3_overlap_re05 | 0.7 | 0.4129 |
| S3_overlap_re05 | 0.8 | 0.4719 |
| S3_overlap_re05 | 0.9 | 0.5308 |
| S4_overlap_re10 | 0.0 | 0.0002 |
| S4_overlap_re10 | 0.1 | 0.0592 |
| S4_overlap_re10 | 0.2 | 0.1182 |
| S4_overlap_re10 | 0.3 | 0.1771 |
| S4_overlap_re10 | 0.4 | 0.2361 |
| S4_overlap_re10 | 0.5 | 0.2951 |
| S4_overlap_re10 | 0.6 | 0.3541 |
| S4_overlap_re10 | 0.7 | 0.4130 |
| S4_overlap_re10 | 0.8 | 0.4720 |
| S4_overlap_re10 | 0.9 | 0.5310 |
| S5_smallK | 0.0 | 0.0000 |
| S5_smallK | 0.1 | 0.0511 |
| S5_smallK | 0.2 | 0.1023 |
| S5_smallK | 0.3 | 0.1534 |
| S5_smallK | 0.4 | 0.2046 |
| S5_smallK | 0.5 | 0.2557 |
| S5_smallK | 0.6 | 0.3069 |
| S5_smallK | 0.7 | 0.3580 |
| S5_smallK | 0.8 | 0.4092 |
| S5_smallK | 0.9 | 0.4603 |
| S6_largeK | 0.0 | 0.0000 |
| S6_largeK | 0.1 | 0.0660 |
| S6_largeK | 0.2 | 0.1321 |
| S6_largeK | 0.3 | 0.1981 |
| S6_largeK | 0.4 | 0.2642 |
| S6_largeK | 0.5 | 0.3302 |
| S6_largeK | 0.6 | 0.3963 |
| S6_largeK | 0.7 | 0.4623 |
| S6_largeK | 0.8 | 0.5284 |
| S6_largeK | 0.9 | 0.5944 |

## What would falsify this

- Simulated E[rho_hat|G] *above* rho (over- rather than under-estimation) at r_e = 0 kills P1 and the paper's framing with it.
- Straight-line residuals above 0.01 kill P3, and the bias would then have to be modelled as genuinely rho-dependent.
- Simulated slope differing from 1/sqrt(a1 a2) by more than ~10% relative means the denominator-inflation account (P2) is incomplete.
