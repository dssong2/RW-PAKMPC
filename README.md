# RW-PAKMPC
PAKMPC but with reaction wheel modeling

## Reaction wheel model

```matlab
build_combined_multibody   % builds pakmpc_rw_combined.slx, runs the free-drift check, ideal case and wheel case
build_rw_array('check')    % rebuilds rw_array_lib.slx and runs the block's four checks
```

| File | What it does |
|---|---|
| `rw_params.m` | All wheel values. `rw_params(preset)` with `'realistic_hr16'` (default, Honeywell HR16-50), `'hr16_100'`, or `'high_torque'` (not a real wheel). Also the wheel geometry, friction, torque error, and the mass and locked inertia the wheels add to the spacecraft. |
| `build_rw_array.m` | Builds the Reaction Wheel Array block in `rw_array_lib.slx`, and with `'check'` runs its checks. |
| `build_combined_multibody.m` | Servicer plus target with the wheel block. Set `rw_preset` at the top. |

The `.slx` files are not tracked. The scripts rebuild them.

Datasheet: Honeywell, "Constellation Series Reaction Wheels", N61-0045-001-001, Dec 2003. HR16: 50/75/100 N m s at 6000 rpm, 0.1 to 0.2 N m nominal torque, unit mass 9/10.4/12 kg.

### Reaction Wheel Array block

Put `rw = rw_params();` in the base workspace and copy the block from `rw_array_lib.slx`.

| Port | Type | Meaning |
|---|---|---|
| `u` | 3x1 in | Body torque command, N m. Sampled every `rw.ts`. |
| `kn` | 1x1 in | Null-space wheel torque, N m. |
| `B` | frame | Connect to the frame the wheels are mounted on. Spin axes are the columns of `rw.wt` in that frame. |
| `p` | 4x1 out | Wheel speeds relative to the body, rad/s. |
| `h` | 4x1 out | Wheel momentum `jw*p`, N m s. |
| `tau` | 3x1 out | Torque on the body, `-W_T*M`, N m. |
| `sat` | 4x1 out | 0 free, 1 torque limit, 2 momentum limit. |
| `pwr` | 4x1 out | Motor power `M.*p`, W. |
| `M` | 4x1 out | Applied motor torque per wheel after the torque error and limits, N m. |

Applied torque, once per sample: `M = -W_T^+ u + kn n - (noise + bias + drift)`, clipped to `taumax`. A wheel being pushed further from zero speed gets at most `(hmax - jw*|p|)/ts`, so it stops at `hmax` instead of passing it. The rotors add `rw.mass` and `rw.ilock` to the spacecraft, so subtract them from the hull.
