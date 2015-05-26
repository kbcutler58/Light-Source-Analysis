function array_in = read_dstBrt(dst_dir, brt_dir, idx_dst, idx_brt)

if exist('idx_brt', 'var')~= 1, idx_brt = 7; end
if exist('idx_dst', 'var')~= 1, idx_dst = 3; end

% read in breast binary mask
pat = dir([brt_dir filesep '*.tif']);
[height, width] = size(imread([brt_dir filesep pat(1).name]));

array_in = zeros(height, width, length(pat));

brtlgh = length(pat);
brtbeg_sl = int8(str2double(pat(1).name(1:3)));

for m = 1 : length(pat),
    array_in(1:height,:,brtlgh - m + 1) = imread([brt_dir filesep ...
        pat(m).name]);
end;

array_in(array_in ~= 0) = idx_brt;

% read in fibro binary mask
pat = dir([dst_dir filesep '*.tif']);
[height, width] = size(imread([dst_dir filesep pat(1).name]));

array_in_dst = zeros(height, width, brtlgh);
dstbeg_sl = int8(str2double(pat(1).name(1:3)));

for m = 1 : length(pat),
    array_in_dst(1:height, :, ...
        brtlgh - abs(brtbeg_sl - dstbeg_sl) + 1 - m) = ...
        imread([dst_dir filesep pat(m).name...
        ]);
end;
% set fibro as 3
array_in(array_in_dst ~= 0) = idx_dst;