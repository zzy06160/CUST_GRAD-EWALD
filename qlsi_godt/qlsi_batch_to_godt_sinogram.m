function sino = qlsi_batch_to_godt_sinogram(cfg, shear_direction)
%QLSI_BATCH_TO_GODT_SINOGRAM  Convert multi-angle QLSI images into GODT sinogram.
%
% sino = qlsi_batch_to_godt_sinogram(cfg, shear_direction)
%
% shear_direction:
%   'x' -> use phiDx and reconstruct approximately dn/dx
%   'y' -> use phiDy and reconstruct approximately dn/dy
%
% The output .mat file is compatible with run_qlsi_godt.m and the GradRytov
% branch of CUST_GRAD-EWALD.

arguments
    cfg struct
    shear_direction char {mustBeMember(shear_direction, {'x','y'})} = 'x'
end

nproj = numel(cfg.input.interferogram_files);
if nproj == 0
    error('cfg.input.interferogram_files is empty.');
end
if numel(cfg.input.reference_files) ~= nproj
    error('reference_files must have the same length as interferogram_files.');
end

[rx, ry] = get_illumination_directions(cfg, nproj);

SINOph = [];
SINOamp = [];
qc = struct();

for j = 1:nproj
    Itf = read_qlsi_image(cfg.input.interferogram_files{j}, cfg.input.mat_variable);
    Ref = read_qlsi_image(cfg.input.reference_files{j}, cfg.input.mat_variable);

    out = qlsi_extract_projection(Itf, Ref, cfg);

    switch shear_direction
        case 'x'
            ph = out.phiDx;
        case 'y'
            ph = out.phiDy;
    end

    ph = real(ph);
    ph = crop_even_square(ph, cfg.output.square_crop_size);
    if cfg.output.subtract_projection_mean
        ph = ph - mean(ph(:), 'omitnan');
    end

    if isempty(SINOph)
        [Ny, Nx] = size(ph);
        if Ny ~= Nx
            error('Projection must be square after cropping.');
        end
        SINOph = zeros(Nx, Nx, nproj, 'single');
        SINOamp = zeros(Nx, Nx, nproj, 'single');
        if cfg.output.save_qc
            qc.OPD = zeros(Nx, Nx, nproj, 'single');
            qc.Tabs = zeros(Nx, Nx, nproj, 'single');
        end
    end

    if ~isequal(size(ph), size(SINOph(:,:,1)))
        error('Projection %d has inconsistent size after cropping.', j);
    end

    SINOph(:,:,j) = single(ph);
    % Transparent weak-scattering first version: amplitude derivative unavailable -> zero.
    SINOamp(:,:,j) = single(zeros(size(ph)));

    if cfg.output.save_qc
        qc.OPD(:,:,j) = single(crop_even_square(out.OPD, cfg.output.square_crop_size));
        qc.Tabs(:,:,j) = single(abs(crop_even_square(out.T, cfg.output.square_crop_size)));
    end
end

lambda_um = cfg.qlsi.lambda_m * 1e6;
NA = cfg.system.NA;

sino_params = zeros(5, nproj);
sino_params(1,:) = rx(:).';
sino_params(2,:) = ry(:).';
sino_params(3,:) = lambda_um .* ones(1, nproj);
sino_params(4,:) = NA .* ones(1, nproj);
sino_params(5,:) = ones(1, nproj); % 1 = transmission

if isempty(cfg.system.dx_um)
    dx = cfg.qlsi.pixel_size_m / cfg.qlsi.relay_zoom / cfg.system.magnification * 1e6;
else
    dx = cfg.system.dx_um;
end

if isempty(cfg.system.shear_um)
    shear = 2 * cfg.qlsi.distance_m * cfg.qlsi.lambda_m / cfg.qlsi.Gamma_m / cfg.system.magnification * 1e6;
else
    shear = cfg.system.shear_um;
end

n_immersion = cfg.system.n_immersion;
geometry = cfg.recon.geometry;
thetay = [];
Fpmask = [];

% Useful scalar aliases for old EWALD conventions.
lambda = lambda_um;
cam_pix = cfg.qlsi.pixel_size_m * 1e6;
downsampling = 1;
M = cfg.system.magnification;

output_folder = cfg.output.folder;
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
output_file = fullfile(output_folder, sprintf('%s_%s.mat', cfg.output.prefix, shear_direction));

save(output_file, 'SINOph', 'SINOamp', 'sino_params', 'shear', 'dx', ...
    'n_immersion', 'geometry', 'thetay', 'Fpmask', 'lambda', 'NA', ...
    'cam_pix', 'downsampling', 'M', 'cfg', 'qc', '-v7.3');

sino = struct();
sino.output_file = output_file;
sino.SINOph = SINOph;
sino.SINOamp = SINOamp;
sino.sino_params = sino_params;
sino.shear = shear;
sino.dx = dx;
sino.n_immersion = n_immersion;
sino.geometry = geometry;
sino.direction = shear_direction;
sino.qc = qc;

fprintf('[QLSI-GODT] Saved %s sinogram: %s\n', shear_direction, output_file);
fprintf('[QLSI-GODT] size = %d x %d x %d, dx = %.6g um, shear = %.6g um\n', ...
    size(SINOph,1), size(SINOph,2), size(SINOph,3), dx, shear);
end

function [rx, ry] = get_illumination_directions(cfg, nproj)
if ~isempty(cfg.illumination.rx) && ~isempty(cfg.illumination.ry)
    rx = cfg.illumination.rx(:).';
    ry = cfg.illumination.ry(:).';
elseif cfg.illumination.use_circular_scan
    if cfg.illumination.Nproj ~= nproj
        error('cfg.illumination.Nproj must match the number of image files.');
    end
    phi = cfg.illumination.phi0_rad + linspace(0, 2*pi, nproj + 1);
    phi(end) = [];
    sin_theta = cfg.illumination.NA_illum / cfg.system.n_immersion;
    rx = sin_theta * cos(phi);
    ry = sin_theta * sin(phi);
else
    error('Provide cfg.illumination.rx/ry or enable cfg.illumination.use_circular_scan.');
end

if numel(rx) ~= nproj || numel(ry) ~= nproj
    error('Illumination direction arrays must match the number of projections.');
end
if any(sqrt(rx.^2 + ry.^2) >= 1)
    error('Invalid illumination direction: sqrt(rx^2+ry^2) must be < 1.');
end
end

function img = read_qlsi_image(file, mat_variable)
[~,~,ext] = fileparts(file);
ext = lower(ext);

switch ext
    case {'.tif','.tiff','.png','.jpg','.jpeg','.bmp'}
        img = imread(file);
        if ndims(img) == 3
            img = rgb2gray(img);
        end
        img = double(img);
    case {'.txt','.csv'}
        img = readmatrix(file);
    case '.mat'
        s = load(file);
        if ~isempty(mat_variable)
            if ~isfield(s, mat_variable)
                error('MAT file %s does not contain variable %s.', file, mat_variable);
            end
            img = s.(mat_variable);
        else
            names = fieldnames(s);
            img = [];
            for k = 1:numel(names)
                candidate = s.(names{k});
                if isnumeric(candidate) && ismatrix(candidate)
                    img = candidate;
                    break;
                end
            end
            if isempty(img)
                error('No numeric 2D matrix found in MAT file %s.', file);
            end
        end
        img = double(img);
    otherwise
        error('Unsupported image extension: %s', ext);
end
end

function out = crop_even_square(img, desired_size)
img = real(img);
[Ny, Nx] = size(img);
N = min(Nx, Ny);
N = 2 * floor(N/2);
if ~isempty(desired_size)
    desired_size = 2 * floor(desired_size/2);
    if desired_size > N
        error('Requested square_crop_size exceeds image size.');
    end
    N = desired_size;
end
cy = floor(Ny/2) + 1;
cx = floor(Nx/2) + 1;
y1 = cy - N/2;
y2 = cy + N/2 - 1;
x1 = cx - N/2;
x2 = cx + N/2 - 1;
out = img(y1:y2, x1:x2);
end
