# QLSI-GODT adapter

This folder contains a minimal adapter layer that combines usable parts of two repositories:

- `CUST_CGMprocess`: QLSI / cross-grating interferogram demodulation.
- `CUST_GRAD-EWALD`: gradient-domain ODT reconstruction with `GradRytov` support.

The intended pipeline is:

```text
QLSI interferograms at multiple illumination angles
        -> qlsi_extract_projection.m
        -> qlsi_batch_to_godt_sinogram.m
        -> Sinogram_QLSI_GODT_x.mat / Sinogram_QLSI_GODT_y.mat
        -> run_qlsi_godt.m
        -> 3D RI-gradient reconstruction
```

## Important distinction

This adapter does **not** use the integrated OPD phase as the main input for tomography. It extracts the shearing phase difference

```text
phi_delta(x,y) = phi(x + Delta, y) - phi(x,y)
```

and sends it to the `GradRytov` branch of the EWALD reconstruction. This is the key difference between a conventional `QLSI -> phase integration -> ODT` workflow and a `QLSI -> GODT` workflow.

## Files

| File | Role |
|---|---|
| `qlsi_godt_config_template.m` | Template for experimental parameters and file lists. |
| `qlsi_extract_projection.m` | Single-projection QLSI demodulation. It extracts shearing phase differences and optional OPD QC data. |
| `qlsi_batch_to_godt_sinogram.m` | Batch conversion from QLSI image files to a GODT-compatible sinogram `.mat` file. |
| `run_qlsi_godt.m` | Non-GUI wrapper for running `CUST_GRAD-EWALD` on the generated sinogram. |
| `qlsi_validate_sinogram.m` | Basic size, unit, NaN, angle and shear checks before reconstruction. |

## Minimal usage

1. Edit `qlsi_godt_config_template.m` and save it as a project-specific config script, for example `config_waveguide_633nm.m`.
2. Fill in image paths, reference paths, illumination directions, grating parameters and microscope magnification.
3. Run:

```matlab
cfg = qlsi_godt_config_template();
sinoX = qlsi_batch_to_godt_sinogram(cfg, 'x');
qlsi_validate_sinogram(sinoX);
[RECONx, infoX] = run_qlsi_godt(sinoX.output_file);
```

For the orthogonal gradient direction:

```matlab
sinoY = qlsi_batch_to_godt_sinogram(cfg, 'y');
qlsi_validate_sinogram(sinoY);
[RECONy, infoY] = run_qlsi_godt(sinoY.output_file);
```

## Required input convention

The generated sinogram contains:

- `SINOph`: shearing phase difference in rad, size `(Nx, Ny, Nproj)`.
- `SINOamp`: amplitude derivative term. The first version sets this to zero for transparent weak-scattering samples.
- `sino_params(1,:)`: dimensionless illumination x-direction component, usually `sin(theta).*cos(phi)`.
- `sino_params(2,:)`: dimensionless illumination y-direction component, usually `sin(theta).*sin(phi)`.
- `sino_params(3,:)`: vacuum wavelength in micrometres.
- `sino_params(4,:)`: detection-objective NA.
- `sino_params(5,:)`: measurement mode, `1` for transmission.
- `shear`: effective shear in the object plane, in micrometres.
- `dx`: object-plane pixel size, in micrometres.
- `n_immersion`: refractive index of the immersion/background medium.
- `geometry`: usually `'fixed'` for stationary sample and scanned illumination.

## Critical checks

Before trusting the result, verify:

1. `SINOph` is an even square stack.
2. `shear` and `dx` are in micrometres.
3. Illumination direction components satisfy `sqrt(rx.^2 + ry.^2) < 1`.
4. The sign of the reconstructed gradient is correct using a known sample or known wavefront tilt.
5. The selected Fourier order and its 90-degree rotated counterpart are stable for all projections.

## Limitations of this first adapter

- Amplitude finite-difference input is set to zero by default. This is acceptable for transparent weak-scattering waveguide samples as a first approximation, but it can cause out-of-focus artifacts for samples with non-negligible amplitude contrast.
- The first-order Fourier crop must be provided or verified. Automatic order detection is intentionally not used in the first version because a wrong order silently corrupts all projections.
- The output is a RI-gradient volume, not an absolute RI volume.
