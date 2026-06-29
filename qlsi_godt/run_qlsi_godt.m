function [RECON, info] = run_qlsi_godt(sinogram_file, overrides)
%RUN_QLSI_GODT  Non-GUI wrapper for GRAD-EWALD GradRytov reconstruction.
%
% [RECON, info] = run_qlsi_godt(sinogram_file)
% [RECON, info] = run_qlsi_godt(sinogram_file, overrides)
%
% This is adapted from GODT.m, but removes the uigetfile dependency and keeps
% all reconstruction parameters explicit. It assumes the sinogram was generated
% by qlsi_batch_to_godt_sinogram.m.

arguments
    sinogram_file char
    overrides struct = struct()
end

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(genpath(repo_root));

if ~exist(sinogram_file, 'file')
    error('Sinogram file not found: %s', sinogram_file);
end

load(sinogram_file); %#ok<LOAD>

% Direct inversion preset from GRAD-EWALD.
DI;

% Reconstruction defaults mirroring GODT.m.
ROI_crop_z = get_override(overrides, 'ROI_crop_z', get_existing_or_default('ROI_crop_z', 2));
limit_resolution_z = get_override(overrides, 'limit_resolution_z', get_existing_or_default('limit_resolution_z', 0.5));
sinogram_subset_factor = get_override(overrides, 'sinogram_subset_factor', get_existing_or_default('sinogram_subset_factor', 1));
projection_crop_factor = get_override(overrides, 'projection_crop_factor', get_existing_or_default('projection_crop_factor', 1));
resample_projections = get_override(overrides, 'resample_projections', get_existing_or_default('resample_projections', false));
plots = get_override(overrides, 'plots', get_existing_or_default('plots', '1'));
save_reconstruction = get_override(overrides, 'save_reconstruction', 0);

% Load missing/legacy defaults used by FDT and Preprocess.
defaults;

% Enforce GODT/gradient mode after defaults, because defaults.m sets Approx='Rytov'.
Approx = 'GradRytov';

if ~exist('SINOamp', 'var') || isempty(SINOamp)
    SINOamp = zeros(size(SINOph), 'single');
    warning('SINOamp missing. Assuming zero amplitude derivative.');
end

if ~exist('shear', 'var') || isempty(shear) || isnan(shear)
    error('Variable shear [um] is required for quantitative QLSI-GODT.');
end

if ~exist('dx', 'var') || isempty(dx) || isnan(dx)
    error('Variable dx [um/pixel] is required for quantitative QLSI-GODT.');
end

if exist('n_immersion', 'var') && ~exist('n_imm', 'var')
    n_imm = n_immersion;
end
if ~exist('n_imm', 'var') || isempty(n_imm)
    error('n_imm or n_immersion is required.');
end

if ~exist('geometry', 'var') || isempty(geometry)
    geometry = 'fixed';
end
if ~exist('thetay', 'var')
    thetay = [];
end
if ~exist('Fpmask', 'var')
    Fpmask = [];
end
if ~exist('do_NNC', 'var')
    do_NNC = 1;
end
if ~exist('n_obj', 'var')
    n_obj = [];
end

[SINOamp_reduced, SINOph_reduced, sino_params, dx, ...
 projection_padding_xy, N_projection_padded, ...
 x_tv, n_obj_res, plots, Fpmask, projection_downsample_factor, N_SINO_y] = ...
    Preprocess(SINOamp, SINOph, geometry, sino_params, resample_projections, dx, projection_crop_factor, Fpmask, ...
               n_obj, thetay, sinogram_subset_factor, n_imm, Masking, N_CP, nCPi, ...
               projection_padding_xy, plots, do_NNC); %#ok<ASGLU>

[RECON, dx_out_xy, dx_out_z, nGPi_done, N_Kspace_xy_padded, KO, KOi, ...
 RMAEtab, RRMSEtab, RMADtab, RRMSDtab] = ...
    FDT(SINOamp_reduced, SINOph_reduced, sino_params, thetay, ...
        n_imm, dx, ...
        shear, ...
        geometry, Approx, interpFp, Ramp, do_NNC, ...
        projection_padding_xy, Kspace_padding, N_projection_padded, Kspace_oversampling_z, ROI_crop_z, limit_resolution_z, ...
        nGPi, epsi, relaxGP, relaxM, ...
        x_tv, n_obj_res, ...
        plots, Fpmask); %#ok<ASGLU>

info = struct();
info.sinogram_file = sinogram_file;
info.Approx = Approx;
info.geometry = geometry;
info.shear_um = shear;
info.input_dx_um = dx;
info.output_dx_xy_um = dx_out_xy;
info.output_dx_z_um = dx_out_z;
info.n_imm = n_imm;
info.nGPi_done = nGPi_done;
info.N_Kspace_xy_padded = N_Kspace_xy_padded;
info.RMAEtab = RMAEtab;
info.RRMSEtab = RRMSEtab;
info.RMADtab = RMADtab;
info.RRMSDtab = RRMSDtab;

if save_reconstruction
    [folder, name] = fileparts(sinogram_file);
    recon_file = fullfile(folder, [name '_RECON.mat']);
    save(recon_file, 'RECON', 'info', '-v7.3');
    info.recon_file = recon_file;
end

if contains(plots, '1')
    figure;
    midz = round(size(RECON,3)/2);
    imagesc(real(RECON(:,:,midz)).');
    axis image;
    colorbar;
    title('QLSI-GODT reconstruction, central XY slice');
    xlabel('x'); ylabel('y');
end
end

function value = get_override(overrides, field, default_value)
if isfield(overrides, field) && ~isempty(overrides.(field))
    value = overrides.(field);
else
    value = default_value;
end
end

function value = get_existing_or_default(varname, default_value)
% This helper intentionally does not inspect caller workspace. It preserves readability
% in the wrapper and keeps defaults local.
value = default_value;
end
