function report = qlsi_validate_sinogram(sino_or_file)
%QLSI_VALIDATE_SINOGRAM  Basic checks for a QLSI-GODT sinogram.
%
% report = qlsi_validate_sinogram(sino_struct)
% report = qlsi_validate_sinogram('Sinogram_QLSI_GODT_x.mat')

if ischar(sino_or_file) || isstring(sino_or_file)
    s = load(char(sino_or_file));
else
    s = sino_or_file;
end

required = {'SINOph','SINOamp','sino_params','shear','dx','n_immersion','geometry'};
for k = 1:numel(required)
    if ~isfield(s, required{k})
        error('Missing required field/variable: %s', required{k});
    end
end

SINOph = s.SINOph;
SINOamp = s.SINOamp;
sino_params = s.sino_params;

[Nx, Ny, nproj] = size(SINOph);
report = struct();
report.Nx = Nx;
report.Ny = Ny;
report.nproj = nproj;
report.is_square = (Nx == Ny);
report.is_even = (mod(Nx,2) == 0) && (mod(Ny,2) == 0);
report.has_nan_phase = any(isnan(SINOph(:)));
report.has_inf_phase = any(isinf(SINOph(:)));
report.has_nan_amp = any(isnan(SINOamp(:)));
report.shear_um = s.shear;
report.dx_um = s.dx;
report.n_immersion = s.n_immersion;
report.geometry = s.geometry;
report.phase_min = min(SINOph(:));
report.phase_max = max(SINOph(:));
report.phase_std = std(double(SINOph(:)), 0, 'omitnan');

if ~isequal(size(SINOamp), size(SINOph))
    error('SINOamp and SINOph must have the same size.');
end
if ~report.is_square
    error('FDT requires square projections. Current size: %d x %d.', Nx, Ny);
end
if ~report.is_even
    error('Use even projection size. Current size: %d x %d.', Nx, Ny);
end
if size(sino_params,1) ~= 5 || size(sino_params,2) ~= nproj
    error('sino_params must be 5 x Nproj.');
end
if report.has_nan_phase || report.has_inf_phase
    error('SINOph contains NaN or Inf.');
end
if isempty(s.shear) || ~isscalar(s.shear) || ~isfinite(s.shear) || s.shear <= 0
    error('shear must be a positive scalar in micrometres.');
end
if isempty(s.dx) || ~isscalar(s.dx) || ~isfinite(s.dx) || s.dx <= 0
    error('dx must be a positive scalar in micrometres per pixel.');
end

rx = sino_params(1,:);
ry = sino_params(2,:);
illum_r = sqrt(rx.^2 + ry.^2);
report.max_illum_direction_norm = max(illum_r);
report.min_illum_direction_norm = min(illum_r);
report.max_illumination_angle_deg = asind(min(0.999999, max(illum_r)));

if any(illum_r >= 1)
    error('Invalid illumination direction: sqrt(rx^2 + ry^2) must be < 1 before GRAD-EWALD scaling.');
end
if any(sino_params(5,:) ~= 1)
    warning('This adapter was designed for transmission mode. Found sino_params(5,:) values other than 1.');
end

fprintf('[QLSI-GODT validation]\n');
fprintf('  size: %d x %d x %d\n', Nx, Ny, nproj);
fprintf('  dx: %.6g um/pixel\n', s.dx);
fprintf('  shear: %.6g um\n', s.shear);
fprintf('  n_immersion: %.6g\n', s.n_immersion);
fprintf('  max illumination direction norm: %.6g\n', report.max_illum_direction_norm);
fprintf('  phase range: [%.6g, %.6g] rad\n', report.phase_min, report.phase_max);
fprintf('  phase std: %.6g rad\n', report.phase_std);
fprintf('  status: basic checks passed\n');
end
