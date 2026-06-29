function out = qlsi_extract_projection(Itf, Ref, cfg)
%QLSI_EXTRACT_PROJECTION  Extract QLSI shearing phase data for GODT.
%
% out = qlsi_extract_projection(Itf, Ref, cfg)
%
% Inputs
%   Itf, Ref : interferogram and reference image, same size.
%   cfg      : configuration struct from qlsi_godt_config_template().
%
% Outputs
%   out.phiDx, out.phiDy : shearing phase differences [rad] in lab x/y axes.
%   out.DWx, out.DWy     : CGM OPD-gradient-like quantities used for QC.
%   out.OPD              : integrated OPD [m], for QC only; do not use as GODT input.
%   out.T                : normalized intensity image.
%   out.crop             : crop metadata used for demodulation.
%
% This function is a self-contained extraction of the usable QLSI demodulation
% logic from CUST_CGMprocess. It intentionally exposes phiDx/phiDy before OPD
% integration so that GODT can operate in the gradient/shearing domain.

Itf = double(Itf);
Ref = double(Ref);

if ~isequal(size(Itf), size(Ref))
    error('Itf and Ref must have the same size.');
end

[Ny, Nx] = size(Itf);
validate_crop(cfg.crop, Nx, Ny);

FItf = fftshift(fft2(Itf));
FRef = fftshift(fft2(Ref));

crop1 = normalise_crop(cfg.crop, Nx, Ny);
crop2 = rotate_crop_90(crop1, Nx, Ny, cfg.crop.rotate90_sign);

H = cell(2,1);
Href = cell(2,1);
[xx, yy] = meshgrid(1:Nx, 1:Ny);

for ii = 1:2
    if ii == 1
        c = crop1;
    else
        c = crop2;
    end
    R2C = (xx - c.x).^2 ./ c.Rx.^2 + (yy - c.y).^2 ./ c.Ry.^2;
    mask = R2C < 1;

    FItfc = FItf .* mask;
    FRefc = FRef .* mask;

    shift_y = -round(c.y) + (Ny/2 + 1);
    shift_x = -round(c.x) + (Nx/2 + 1);
    H{ii} = circshift(FItfc, [shift_y, shift_x]);
    Href{ii} = circshift(FRefc, [shift_y, shift_x]);
end

Ix    = ifft2(ifftshift(H{1}));
Iy    = ifft2(ifftshift(H{2}));
Irefx = ifft2(ifftshift(Href{1}));
Irefy = ifft2(ifftshift(Href{2}));

phase_sign = cfg.crop.phase_sign;
phiD1 = phase_sign .* angle(Ix .* conj(Irefx));
phiD2 = phase_sign .* angle(Iy .* conj(Irefy));

ct = cos(crop1.theta_rad);
st = sin(crop1.theta_rad);

% Rotate from grating/order coordinates to lab x/y coordinates.
phiDx = ct .* phiD1 - st .* phiD2;
phiDy = st .* phiD1 + ct .* phiD2;

% CGM OPD-gradient-like outputs for QC. These are not the preferred GODT input.
alpha = cfg.qlsi.Gamma_m / (4*pi*cfg.qlsi.distance_m);
DWx = alpha .* phiDx;
DWy = alpha .* phiDy;

% Optional OPD integration for QC only.
OPD = integrate_gradients_fft(DWx, DWy, cfg.qlsi.pixel_size_m, cfg.qlsi.relay_zoom);

% Normalized intensity from the zero-order/central band.
R0x = crop1.Rx;
R0y = crop1.Ry;
R2C0 = (xx - (Nx/2 + 1)).^2 ./ R0x.^2 + (yy - (Ny/2 + 1)).^2 ./ R0y.^2;
mask0 = R2C0 < 1;
HT = FItf .* mask0;
HTref = FRef .* mask0;
T = ifft2(ifftshift(HT)) ./ ifft2(ifftshift(HTref));

out = struct();
out.phiDx = real(phiDx);
out.phiDy = real(phiDy);
out.phiD1 = real(phiD1);
out.phiD2 = real(phiD2);
out.DWx = real(DWx);
out.DWy = real(DWy);
out.OPD = real(OPD);
out.T = T;
out.crop = crop1;
out.crop2 = crop2;
out.meta.Gamma_m = cfg.qlsi.Gamma_m;
out.meta.distance_m = cfg.qlsi.distance_m;
out.meta.pixel_size_m = cfg.qlsi.pixel_size_m;
out.meta.relay_zoom = cfg.qlsi.relay_zoom;
out.meta.lambda_m = cfg.qlsi.lambda_m;
end

function validate_crop(crop, Nx, Ny)
if ~isfield(crop,'x') || isempty(crop.x) || ~isfield(crop,'y') || isempty(crop.y)
    error('cfg.crop.x and cfg.crop.y must be provided. Select a first Fourier order manually first.');
end
if (~isfield(crop,'R') || isempty(crop.R)) && ...
   ((~isfield(crop,'Rx') || isempty(crop.Rx)) || (~isfield(crop,'Ry') || isempty(crop.Ry)))
    error('Provide cfg.crop.R or both cfg.crop.Rx and cfg.crop.Ry.');
end
if crop.x < 1 || crop.x > Nx || crop.y < 1 || crop.y > Ny
    error('Crop centre is outside the image.');
end
end

function crop = normalise_crop(crop, Nx, Ny)
if ~isfield(crop,'Rx') || isempty(crop.Rx)
    crop.Rx = crop.R;
end
if ~isfield(crop,'Ry') || isempty(crop.Ry)
    crop.Ry = crop.R;
end
if ~isfield(crop,'theta_rad') || isempty(crop.theta_rad)
    cx = Nx/2 + 1;
    cy = Ny/2 + 1;
    crop.theta_rad = atan2(crop.y - cy, crop.x - cx);
end
if ~isfield(crop,'phase_sign') || isempty(crop.phase_sign)
    crop.phase_sign = 1;
end
if ~isfield(crop,'rotate90_sign') || isempty(crop.rotate90_sign)
    crop.rotate90_sign = 1;
end
end

function crop2 = rotate_crop_90(crop1, Nx, Ny, rotate90_sign)
% Rotate the selected Fourier-order crop around the FFT centre.
cx = Nx/2 + 1;
cy = Ny/2 + 1;
dx = crop1.x - cx;
dy = crop1.y - cy;

if rotate90_sign >= 0
    dx2 = -dy;
    dy2 =  dx;
else
    dx2 =  dy;
    dy2 = -dx;
end

crop2 = crop1;
crop2.x = cx + dx2;
crop2.y = cy + dy2;
crop2.theta_rad = crop1.theta_rad + rotate90_sign*pi/2;
end

function OPD = integrate_gradients_fft(DWx, DWy, pixel_size_m, relay_zoom)
[Ny, Nx] = size(DWx);
[kx, ky] = meshgrid(1:Nx, 1:Ny);
kx = kx - Nx/2 - 1;
ky = ky - Ny/2 - 1;
zero_id = (kx == 0) & (ky == 0);
kx(zero_id) = Inf;
ky(zero_id) = Inf;

W0 = ifft2(ifftshift((fftshift(fft2(DWx)) + 1i*fftshift(fft2(DWy))) ./ ...
    (1i*2*pi*(kx/Nx + 1i*ky/Ny))));
OPD = pixel_size_m/relay_zoom * real(W0);
end
