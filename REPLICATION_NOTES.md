# PAKMPC replication notes

MATLAB R2026b replication of Sun et al., "Physics-based Online Adaptive Koopman Model Predictive
Attitude Control for Combined Spacecraft with Dynamic Uncertainties" (arXiv 2609.02534,
IEEE TASE, DOI 10.1109/TASE.2026.3730829), on branch `dev/pakmpc`.

These notes list every choice the paper leaves open (tagged `[guess]` in the code) and every
place where the results differ from the paper. Runtimes are machine-dependent throughout, so only
the ordering between methods is compared.

## What was built

| Step | Paper result | Script | Library |
|---|---|---|---|
| 1 | Fig. 2, nominal model prediction | `fig2_nominal_models.m` | `quat_koopman`, `dcm_koopman` |
| 2 | Fig. 3, Table II, Fig. 4 | `fig3_table2_nominal_mpc.m`, `fig4_monte_carlo.m` | `koopman_mpc`, `rbf_edmd`, `terminal_weight`, `simulate_attitude`, `settling_time` |
| 3 | Figs. 6-11, Table III (PAKMPC) | `fig6_11_table3_adaptive.m` | `koopman_mpc` (`adapt`) |
| 4 | KMAE-AKMPC and INMPC baselines | same script | `kmae_akmpc`, `incremental_nmpc` |
| 5 | Figs. 12-13, robustness sweep | `fig12_13_robustness.m` | `combined_inertia` |
| 6 | Fig. 5, Simscape plant | `fig5_simscape.m` | `simscape_step`, `draw_combined` |

Development used an `ode45` rigid-body plant (`simulate_attitude`). Step 6 confirms that this plant
agrees with Elijah's Simscape Multibody model (`build_pakmpc_multibody.m`): without noise,
the closed-loop attitude matches to 6e-12 and the torque to 1e-9 N m. So every ode45 result
here also holds on the multibody plant.

## Choices the paper leaves open (`[guess]`)

### Models

- **DCM model size** (`dcm_koopman`). The text says 9(N_phi + 1) states while Table II gives 9 N_phi.
  Reference [18] isn't reproduced, so w is kept in the state as in the quaternion model,
  giving 9 N_phi + 12.
- **RBF-EDMD training data** (`rbf_edmd`). Not given in the paper. Used here:
  - 200 trajectories of 100 samples, starting from random attitudes with |w_i| <= 0.02 rad/s.
  - Inputs uniform in [-umax, umax], held over each sample.
  - 50 thin-plate RBF centers drawn from the data; w scaled by 1/0.05.
- **Fig. 2 test conditions** (`fig2_nominal_models`):
  - J_c0 = J_c, with no inertia uncertainty, as in Sec. V-B.
  - w0 = 0, to match Fig. 2's flat start (the Sec. V-A w0 gives error from t = 0).
  - No disturbance; the `dist_on` switch adds it.
  - B is re-evaluated at the model's own predicted state, so the run is a pure open-loop prediction.

### MPC (`koopman_mpc`, `terminal_weight`)

- **Terminal weight P.** This replaces the LMI (42)-(43), which has no solution for this model.
  - Every eigenvalue of the lifted A is 1, and rank [A - I, B0] < n, so (A, B0) isn't stabilizable.
  - Instead, P6 solves the DARE of the linearization of [w; q_v] at the target, and is mapped
    to the lifted state as P = C' T6' P6 T6 C.
  - The paper's κ, H_x, h_x and H_u aren't given, so σ = 208.68 can't be reproduced or used.
  - Without a real P, the 1.5 s horizon barely acts: peak torque is 4.5 N m instead of saturating.
- **Cost.** The cost penalizes z - Phi(x_e), so q -> [1 0 0 0]. Eq. (29a) as printed would drive q0 to 0.
- **Input bound.** |τ_i| <= u_max is mapped to bounds on u = τ - w^x J_c0 w, using the
  current measured w held over the horizon.
- **Nominal B.** B is rebuilt from the measured (w, q) every sample. The `freeze_b` option holds
  B(q0, w0) instead, per Algorithm 1 line 2 and Eq. 38; it made no difference.
- **Terminal set.** The terminal set (28g) is left out of the QP, as in the paper's QP form (34).
- **DCM weights.** The DCM model uses Q = 1e2 on all 12 outputs.
- **Λ and γ_e.** Λ = diag(λ^N_w, ..., λ) follows the sum in (25); the printed Λ has N_w - 1 entries.
  The "λμ²" term in the QP is read as γ_e μ².
- **Adam update.** One Adam step per sample. u_(k-1) is the applied (saturated) torque, mapped back
  through Eq. (8).

### Plant and scenario (`pakmpc_params`, `simulate_attitude`)

- **Measurement noise.** The paper writes N(0, 1e-4 I3). This is read as standard deviation 1e-4
  (`p.noise_var = 1e-8`).
  - Read as a variance (std 0.01), the Fig. 8 torque is buried in noise and ||ΔB|| drifts up
    instead of decaying as in Fig. 11.
  - Elijah's `build_pakmpc_multibody.m` uses variance 1e-4; `fig5_simscape` overrides it.
- **Target torque reference.** Eq. (59) "q_vt ... of the initial quaternion of the target" is
  ambiguous. The adaptive scripts regulate the target toward the identity attitude
  (`target_ref = 'identity'`), which reproduces Table III's steady-state ||q_v||.
  - Holding the initial attitude (`'initial'`, as in Elijah's builder) leaves about 7e-3.
  - `fig5_simscape` overrides the builder's target block to match.
- **Steady-state window.** Steady-state errors are averaged over t = 190-200 s, the window of the
  Fig. 6 insets.
- **Settling time** (`settling_time`). The paper says "within a 5% error margin" without saying
  of what; here it is 5% of the initial ||q_v||.
- **Fig. 4 seed and quaternion sign.** Random seed 42, and initial quaternions with q0 >= 0
  (shortest rotation).

### Baselines (step 4)

The paper describes each baseline in one sentence, so their tuning is entirely chosen here.

- **KMAE-AKMPC** (`kmae_akmpc`, after [31]):
  - Starts from the RBF-EDMD model trained on the nominal inertia J_c0 = J_s.
  - Updates the full [A B] by RLS, continuing the offline fit (P0 = inv(Gram)) with forgetting factor 0.995.
  - Runs the same linear MPC as PAKMPC, without the interpolated initial state.
- **INMPC** (`incremental_nmpc`, after [41]):
  - INDI-type incremental model: the measured angular acceleration, low-pass filtered (α = 0.5),
    plus J_c0^-1 times the change in torque, with full quaternion kinematics.
  - Same cost and terminal weight as the other controllers.
  - Solved with fmincon SQP (30 iterations), warm-started from the previous solution.

### Robustness sweep and Simscape (steps 5-6)

- **Combined inertia** (`combined_inertia`). The paper gives J_c but not how it is assembled.
  - Rebuilt from the servicer, the target and the docking geometry of `build_pakmpc_multibody`.
  - The remaining inertia (arms, panels, ring) is held fixed while the target's mass and
    inertia are scaled.
  - Returns the paper's J_c exactly at 100%.
- **Target gains.** K_pt = 0.02 J_t and K_dt = 0.05 J_t follow the scaled J_t, per Eq. (59).
- **Noise seed.** Every sweep run uses the same noise seed.
- **Fig. 5 snapshots** (`draw_combined`).
  - The paper's are Mechanics Explorer screenshots, which can't be captured headless, so the
    snapshots are drawn with MATLAB patches instead.
  - They use Elijah's geometry, itself a guess: his target box is slightly larger than the
    servicer, while the paper's servicer is much larger.

## Mismatches with the paper

### Main result: no gap between PAKMPC and non-adaptive MPC (Figs. 6-11, Table III)

The paper's non-adaptive MPC overshoots heavily and keeps saturating, at about 3.6x PAKMPC's
control effort. Here the two are nearly identical:

| | Steady-state ‖q_v‖ [paper Table III] | Steady-state ‖w‖ [Table III] | Effort ∫‖τ‖² dt [paper Fig. 12b, nominal cell] |
|---|---|---|---|
| PAKMPC | 1.29e-3 [2.1e-3] | 4.95e-4 [2.2e-4] | 1.00e4 [about 0.45e4] |
| Non-adaptive | 1.29e-3 [2.8e-3] | 4.95e-4 [2.1e-4] | 1.05e4 [about 1.8e4] |

Why: with ρ = 1e-7, Adam moves ΔB by at most ρ per step per entry.
- ||ΔB|| therefore stays near 1e-5 to 4e-5, while correcting the 2.5x inertia error would need about 1.8e-4.
- The closed loop is dominated by the fixed terminal weight over the short 1.5 s horizon.
  Even an oracle ΔB, which makes the prediction exact, barely changes the response.

None of these readings recover the gap:
- ΔB accumulated into B every sample (Algorithm 1 line 11): the loop diverges.
- B frozen at B(q0, w0): no change.
- ρ raised from 1e-6 to 1e-4: worse, as ΔB fits noise; at 1e-4 the loop doesn't settle.
- Noise as a variance: no gap.
- Target holding its initial attitude: no gap.
- Weaker non-adaptive baselines.

The arXiv TeX source has no hidden parameters and no LMI setup. Its unused draft figure
`heatmaps.pdf` is an earlier PAKMPC vs non-adaptive comparison, and it shows similar steady-state
errors for the two (about 2.5e-3 vs 3e-3), much like these results. The published gap most likely
depends on the undocumented terminal weight and baseline settings.

### Other differences

| Result | Paper | Here |
|---|---|---|
| Fig. 2 error at 300 s, quaternion N_phi = 5-9 | about 0.4 | 0.02-0.03 (adding τ_d to the plant gives 0.5, but loses the flat start) |
| Fig. 2 error at 300 s, DCM N_phi = 7-9 | about 0.4 | about 0.1 |
| Fig. 3 quaternion settling | about 40-50 s | about 25 s |
| Fig. 3 DCM | "inverted" torque | not seen |
| Fig. 3 RBF-EDMD | slow, oscillatory | less oscillatory (training data is a guess) |
| Fig. 4 median settling time | about 115 s EDMD, about 30 s analytic models | 27.7 s EDMD, about 18 s analytic models (same ordering and spread pattern) |
| Table II runtime | 16.85 s per 200 s run (quaternion) | 1.35 s (same ordering: quaternion < EDMD ≈ DCM, DCM slows most with N_phi) |
| Fig. 10 μ_k | 1 until about 25 s, then decays | stays below 0.3 |
| Fig. 11 ‖ΔB‖ peak | about 1.5e-5, gone by 10-15 s | about 4e-5, decays over 200 s |
| Table III INMPC ‖q_v‖ | 3.9e-3 (worst) | 4.35e-4 (best): the incremental model measures and cancels the disturbance and target torque directly |
| Table III KMAE-AKMPC ‖q_v‖ | 3.4e-3, worse than PAKMPC | 1.59e-3, worse than PAKMPC (ordering matches) |
| Table III step time | PAKMPC 0.010 < KMAE 0.017 < INMPC 0.035 s | 0.0007 < 0.0009 < 0.0076 s (ordering matches) |
| Fig. 12a | error grows with target inertia, most for KMAE-AKMPC and INMPC (up to 8e-3) | error depends almost only on disturbance scale. Scaling a 75 kg target against a 1080 kg servicer moves J_c by only about ±13% |
| Fig. 12b | non-adaptive about 1.8e4, others about 0.45e4 | all about 1e4, increasing with inertia |
| Fig. 13 | PAKMPC nearest the lower left | PAKMPC and non-adaptive overlap; INMPC most accurate but slowest |

### What does match

- Model structure: the lifted derivative matches the equations to 1e-10 for both analytic models,
  with only the truncated last block differing.
- Fig. 2 trends: higher N_phi gives lower error, and the quaternion model beats DCM at equal N_phi.
- Fig. 3 / Table II / Fig. 4 rankings: the quaternion model is fastest, smallest and most
  consistent; EDMD is slowest with the widest spread.
- Figs. 6-8 shape: torque saturates at the start, attitude settles by about 50 s, and steady-state errors are of the same order.
- Fig. 11 shape: ||ΔB|| rises and then decays, of the same order of magnitude.
- KMAE-AKMPC trails PAKMPC, and runtime orders in Table III and Fig. 13 agree.
- PAKMPC uses slightly less control effort than non-adaptive MPC across the whole sweep.

## Reproducing

All scripts run headless with `matlab -batch "run('<script>.m')"` and write to `results/`.
`fig4_monte_carlo` and `fig12_13_robustness` use `parfor`. `fig5_simscape` needs Simscape
Multibody and builds `pakmpc_multibody.slx` (git-ignored) on first run.
