function [h_out, meshBrt] = disp_dstBrt3D(array_in, option)

if exist('option', 'var')~= 1
    option.dst.idx = 3;
    option.dst.color = [0.95 0.95 0.95];
    option.dst.alpha = .3;
    option.dst.iso = .7;
    option.brt.idx = 7;
    option.brt.color = [0.95 0.95 0.95];
    option.brt.alpha = .2;
    option.brt.iso = .5;
end

if ~isfield(option, 'dst'), option.dst = []; end
if ~isfield(option, 'brt'), option.brt = []; end
if ~isfield(option.dst, 'idx'), option.dst.idx = 3; end
if ~isfield(option.dst, 'color'), option.dst.color = [0.95 0.95 0.95]; end
if ~isfield(option.dst, 'alpha'), option.dst.alpha = .3; end
if ~isfield(option.dst, 'iso'), option.dst.iso = .7; end
if ~isfield(option.brt, 'idx'), option.brt.idx = 7; end
if ~isfield(option.brt, 'color'), option.brt.color = [0.95 0.95 0.95]; end
if ~isfield(option.brt, 'alpha'), option.brt.alpha = .2; end
if ~isfield(option.brt, 'iso'), option.brt.iso = .5; end

pad = 5;
smoothIdxBrt = 13;
smoothIdxDst = 9;
isoValueDst = option.dst.iso;
isoValueBrt = option.brt.iso;

[y x z] = ind2sub(size(array_in), find(array_in));
sy = max(y) - min(y) + 1 + pad *2;
sx = max(x) - min(x) + 1 + pad *2;
sz = max(z) - min(z) + 1 + pad *2;
imgVol = zeros([sy sx sz]);
imgVol(pad + 1: sy - pad, pad + 1: sx - pad, pad + 1: sz - pad) = ...
    array_in(min(y):max(y), min(x): max(x), min(z): max(z));
array_in = imgVol;

if ~isfield(option, 'thick'), thick = 1/size(imgVol, 3); 
else thick = option.thick;
end
if ~isfield(option, 'pixel_xy'), pixel_xy = 1/size(imgVol, 2); 
else pixel_xy = option.pixel_xy;
end

array_in_dst = zeros(size(array_in));
array_in_dst(array_in == option.dst.idx) = 1;
array_in_brt = zeros(size(array_in));
array_in_brt(array_in == option.brt.idx) = 1;

brtdata = smooth3(array_in_brt,'box',smoothIdxBrt);

[x_cent, ~, z_cent] = centroid3D(array_in_brt);
meshBrt = isosurface(brtdata, isoValueBrt);
meshBrt.vertices =  bsxfun(@minus, meshBrt.vertices, ...
    [x_cent max(y) z_cent]);
h.brt = patch(meshBrt,'FaceColor', option.brt.color,...
    'EdgeColor','none');
alpha(h.brt, option.brt.alpha);
% isonormals(brtdata, brtph);

if isfield(option.brt, 'extra') && ~isempty(option.brt.extra)
    h.extra = patch(meshBrt,'FaceColor', ...
        option.brt.extra,...
        'EdgeColor','none');
    alpha(h.extra, 0.15);
%     isonormals(brtdata, brtph_extra); 
end

view(3); axis vis3d tight; 
camlight('local'); 
lighting phong;
material dull
xlabel('x-axis'),ylabel('y-axis'),zlabel('z-axia');

dstdata = smooth3(array_in_dst,'box',smoothIdxDst);
meshDst = isosurface(dstdata, isoValueDst);
meshDst.vertices =  bsxfun(@minus, meshDst.vertices, ...
    [x_cent max(y) z_cent]);
h.dst = patch(meshDst,'FaceColor', option.dst.color,...
    'EdgeColor','none');
alpha(h.dst, option.dst.alpha)
daspect(1./[pixel_xy pixel_xy thick]);

if isfield(option.brt, 'extra') && ~isempty(option.brt.extra)
    h_out = [h.brt h.dst h.extra];
else
    h_out = [h.brt h.dst];
end
