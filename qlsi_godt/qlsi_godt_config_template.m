function cfg = qlsi_godt_config_template()
%QLSI_GODT_CONFIG_TEMPLATE  Configuration template for QLSI-GODT conversion.
%
% Copy this file to a project-specific config, then edit the values.
% Units:
%   - QLSI optical geometry below is entered in metres.
%   - GODT sinogram variables are saved in micrometres where required.
%
% The illumination direction convention follows the GRAD-EWALD preprocessing:
%   sino_params(1:2,:) are dimensionless direction components before scaling,
%   i.e. rx = sin(theta)*cos(phi), ry = sin(theta)*sin(phi).

%% Input images
% Use full paths or paths relative to the MATLAB current folder.
% One interferogram and one reference image are required per illumination angle.
cfg.input.interferogram_files = {
    % 'data/proj_0001.tif'
    % 'data/proj_0002.tif'
};

cfg.input.reference_files = {
    % 'data/ref_0001.tif'
    % 'data/ref_0002.tif'
};

% Supported by read_qlsi_image() inside qlsi_batch_to_godt_sinogram:
% .tif/.tiff/.png/.jpg/.jpeg/.bmp, .txt/.csv, .mat.
% For .mat files, set cfg.input.mat_variable if needed.
cfg.input.mat_variable = '';

%% QLSI / CGM parameters, image-side units
cfg.qlsi.Gamma_m      = 39e-6;   % cross-grating period [m]
cfg.qlsi.distance_m   = 1.0e-3;  % effective grating-to-camera/image distance [m]
cfg.qlsi.pixel_size_m = 6.5e-6;  % physical camera pixel size [m]
cfg.qlsi.relay_zoom   = 1.0;     % relay-lens magnification between grating and sensor description
cfg.qlsi.lambda_m     = 633e-9;  % vacuum wavelength [m]

% Microscope magnification from sample to camera. If the effective sample-plane
% pixel size is already known, set cfg.system.dx_um directly below.
cfg.system.magnification = 50;

% Detection objective NA and immersion/background refractive index.
cfg.system.NA = 0.75;
cfg.system.n_immersion = 1.0;    % use glass/background RI for embedded waveguides when appropriate

% Object-plane sampling [um/pixel]. Leave [] to compute from pixel_size/relay_zoom/magnification.
cfg.system.dx_um = [];

% Effective shear in object plane [um]. Leave [] to use nominal 2*d*lambda/Gamma/M.
% For quantitative GODT, prefer experimental shear calibration.
cfg.system.shear_um = [];

%% Fourier first-order crop
% This first version expects the user to provide a stable first-order crop.
% x,y are in MATLAB image coordinates: x = column index, y = row index.
% R is the circular crop radius in pixels. Rx/Ry can be used for elliptical crop.
% theta_rad is the angle between the selected Fourier order axis and the lab x-axis.
cfg.crop.x = [];          % e.g. 1190
cfg.crop.y = [];          % e.g. 850
cfg.crop.R = [];          % e.g. 70
cfg.crop.Rx = [];
cfg.crop.Ry = [];
cfg.crop.theta_rad = [];  % if empty, computed from selected order position relative to FFT centre

% Phase sign convention. Change to -1 if validation with a known wavefront shows inverted gradients.
cfg.crop.phase_sign = 1;

% Rotation direction for the second order. +1 = counter-clockwise 90 deg; -1 = clockwise 90 deg.
cfg.crop.rotate90_sign = 1;

%% Illumination directions
% Fill one value per projection. These are dimensionless components of illumination direction.
% For a circular scan:
%   theta = asin(NA_illum / n_immersion);
%   rx = sin(theta) * cos(phi);
%   ry = sin(theta) * sin(phi);
cfg.illumination.rx = []; % row vector, length Nproj
cfg.illumination.ry = []; % row vector, length Nproj

% Optional helper for circular scan generation. If rx/ry above are empty and this block is enabled,
% qlsi_batch_to_godt_sinogram will generate rx/ry from it.
cfg.illumination.use_circular_scan = false;
cfg.illumination.Nproj = 61;
cfg.illumination.NA_illum = 0.55;
cfg.illumination.phi0_rad = 0;

%% Sinogram output
cfg.output.folder = 'qlsi_godt/output';
cfg.output.prefix = 'Sinogram_QLSI_GODT';

% Projection preprocessing before saving. Use an even square crop for FDT.
cfg.output.square_crop_size = [];  % e.g. 1024. Leave [] to use largest even centered square.
cfg.output.subtract_projection_mean = true; % removes constant piston from shearing phase difference
cfg.output.save_qc = true;

%% GRAD-EWALD reconstruction parameters
cfg.recon.geometry = 'fixed';
cfg.recon.plots = '1';
cfg.recon.ROI_crop_z = 2;
cfg.recon.limit_resolution_z = 0.5;
cfg.recon.sinogram_subset_factor = 1;
cfg.recon.projection_crop_factor = 1;
cfg.recon.resample_projections = false;

end
