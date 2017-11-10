function rgb =  munsell2rgb(hue, value, chroma, varargin)
% MUNSELL2RGB Convert Munsell color chip coordinates to RGB values
%   MUNSELL2RGB(H, V, C) uses the specified H, V, C coordinates to lookup
%   x, y, Y coordinates which are transformed to XYZ coordinates to sRGB
%   and finally to RGB coordinates using transformation matrices available
%   at Lindbloom.com
%
% Optional inputs: 
%   is_sRGB = true|{false}
%   gamma   = {2.2}
%
%
% Notes: 
% To emulate the [R] code at http://casoilresource.lawr.ucdavis.edu/drupal/node/201, set is_sRGB = true
%
% To emulate the WallKillColor program, set is_sRGB = false, gamma = 2.2. 
% This is the default setting
global is_sRGB

optargs = {false, 2.2}; 

% skip any new inputs if they are empty, and put these defaults into the 
% valuesToUse cell array, and overwrite the ones specified in varargin.
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);

% Place optional args in memorable variable names
[is_sRGB, gamma] = optargs{:};

%% Read Munsell table
% Munsell chroma, CIE x, y, and Y. 
% The chromaticity coordinates were calculated using illuminant C and the CIE 1931 2 degree observer
% http://www.cis.rit.edu/mcsl/online/munsell.php

fid = fopen('all.dat');
muntab = textscan(fid, '%s%f%f%f%f%f',  'Headerlines', 1);
fclose(fid);

%% Put data from table into variables
% Munsell coordinates
H = muntab{1}; % hue
V = muntab{2}; % value/brightness
C = muntab{3}; % chroma/saturation

% xyY coordinates
xyY = [muntab{4:6}];

%% Check input
hueset = unique(H); valset = unique(V); chrset = unique(C); 
if sum(strcmp(hueset, hue)) == 0; error('Unknown hue value specified.'); end 
if value < min(valset)  || value  > max(valset); error('Value out of range.');  end
if chroma < min(chrset) || chroma > max(chrset); error('Chroma out of range.'); end

%% Locate input in table
if any(strcmp(H, hue) & V == value & C == chroma) % If there is an exact match, retrieve it...
    cset = xyY(strcmp(H, hue) & V == value & C == chroma, :);
else   % locate nearby points and interpolate
%     vL = valset(find(valset <= value, 1, 'last')); 
%     vU = valset(find(valset >= value, 1, 'first')); 
%     
%     cL = chrset(find(chrset <= chroma, 1, 'last')); 
%     cU = chrset(find(chrset >= chroma, 1, 'first'));
%     
%     csetL = xyY(strcmp(H, hue) & V == vL & C == cL, :);
%     csetU = xyY(strcmp(H, hue) & V == vU & C == cU, :);
%     
%     % Do simple interpolation
%     cset  = csetL + (csetU - csetL)/2;
    cset(1) = interpolate(H, V, C, xyY, value, chroma, hue, 'x');
    cset(2) = interpolate(H, V, C, xyY, value, chroma, hue, 'y'); 
    cset(3) = interpolate(H, V, C, xyY, value, chroma, hue, 'Y'); 
end
x = cset(1); y = cset(2); Y = cset(3);

%% Convert xyY to XYZ
% http://www.brucelindbloom.com/index.html?Eqn_xyY_to_XYZ.html
% x and y are approx (0,1)
% Y is approx (0,100)

% need manually rescale Y to (0,1)
Y = Y/100.0;

% do the conversion
X = (x .* Y ) ./ y;
% Y = Y;
Z = ( (1 - x - y) .* Y )  ./ y;

% combine to form matrix for simple manipulation
mun_XYZ_C = [X,Y,Z];

% if y == 0, X,Y,Z should then be set to 0
mun_XYZ_C(y == 0, :) = 0;

%% Perform chromatic adaptation 
% from http://www.brucelindbloom.com/index.html?Eqn_ChromAdapt.html)
% NOTE: Lindbloom gives the transpose matrix, but the order of the matrix
% operations is reversed so the matrix when using Lindbloom's matrices,
% transpose them first.
% M_adapt_C_to_D65 = [0.990448, -0.012371, -0.003564; -0.007168, 1.015594, 0.006770; -0.011615, -0.002928, 0.918157];
M_adapt_C_to_D65 = [0.9904476 -0.0071683 -0.0116156; -0.0123712  1.0155950 -0.0029282; -0.0035635  0.0067697  0.9181569]';

% perform the chromatic adaption: convert from C -> D65 using Bradford method
mun_XYZ_D65 = mun_XYZ_C * M_adapt_C_to_D65;


% how different are the two?
difference = (mun_XYZ_D65 - mun_XYZ_C);

%% Convert XYZ (D65) to sRBG (D65), step 1
% http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
% assumes that XYZ is scaled to (0, 1)

% first get the reference primaries transformation matrix from above

% sRGB profile transformation:
M_XYZ_to_sRGB_D65 = [3.24071, -0.969258, 0.0556352; -1.53726, 1.87599, -0.203996; -0.498571, 0.0415557, 1.05707];

% apply the conversion matrix
mun_sRGB_D65 = mun_XYZ_D65 * M_XYZ_to_sRGB_D65;

%% Convert XYZ (D65) to sRGB (D65), step 2 (sRGB, gamma = 2.4) 
% http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_RGB.html
% NOTE: The values output from this code are closer to WallKillColor gamma = 2.2

% the specific function is contingent on the absolute value of r,g,b components
R = mun_sRGB_D65(:,1);
R(R  > 0.0031308) = fun1(R(R  > 0.0031308), gamma);
R(R <= 0.0031308) = fun2(R(R <= 0.0031308), gamma);

G = mun_sRGB_D65(:,2);
G(G >  0.0031308) = fun1(G(G  > 0.0031308), gamma);
G(G <= 0.0031308) = fun2(G(G <= 0.0031308), gamma);

B = mun_sRGB_D65(:,3);
B(B >  0.0031308) = fun1(B(B >  0.0031308), gamma);
B(B <= 0.0031308) = fun2(B(B <= 0.0031308), gamma);

% clip values to range {0,1}
R(R < 0) = 0; R(R > 1) = 1;
G(G < 0) = 0; G(G > 1) = 1;
B(B < 0) = 0; B(B > 1) = 1;
 
rgb = round([R G B] * 255);

%% transformation functions: these are applied on a conditional basis:
function out = fun1(col_comp, gamma) 
global is_sRGB

switch is_sRGB
    case true % If colors are sRGB
        out = 1.055 * ( col_comp .^ ( 1 / 2.4 ) ) - 0.055;
    case false  % ...if colors are not sRGB
        out = col_comp .^ (1/gamma);
end

function out = fun2(col_comp, gamma) 
global is_sRGB

switch is_sRGB
    case true % If colors are sRGB
        out = 12.92 * col_comp;
    case false  % ...if colors are not sRGB
        out = col_comp .^ (1/gamma);
end

function ival = interpolate(H, V, C, xyY, value, chroma, hue, val)
    valset = unique(V); chrset = unique(C);
    subset = sortrows([V(strcmp(hue, H)), C(strcmp(hue, H)), xyY(strcmp(hue, H), strmatch(val, {'x', 'y', 'Y'}))], [1 2]);
    [vm, cm] = meshgrid(valset, chrset); 
    zm = nan(numel(chrset), numel(valset));
    for i = 1:size(subset, 1)
        zm(vm == subset(i,1) & cm == subset(i, 2)) = subset(i, 3);
    end
    ival = interp2(vm, cm, zm, value, chroma); 