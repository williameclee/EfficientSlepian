%% POLARGRIDMASK - Creates grid for the polar cap and domain mask
%
% Syntax
%   [lon, lat, res, weight, mask] = polarGridMask(radius, pLonlat)
%   [lon, lat, res, weight, mask] = polarGridMask(radius, pLonlat, L)
%   [lon, lat, res, weight, mask] = polarGridMask(__, "Name", value)
%
% Input arguments
%   radius - Scalar radius of the polar cap, in degrees
%   pLonlat - Nx2 array of the rotated domain boundary
%       [longitude, latitude] in degrees
%   L (optional) - Bandwidth (maximum angular degree), used to determine
%       grid resolution
%       The default degree is 18.
%   resFactor (name-value) - Resolution scaling factor; higher values
%       produce finer grids
%       The default value is 8
%
% Output arguments
%   lon - Grid longitudes, in degrees
%       Size: [1 x M]
%   lat - Grid latitudes, in degrees
%       Size: [N x 1]
%   res - Grid angular resolution, in degrees
%   weight - Integration weights (solid-angle element)
%       Size: [N x M]
%   mask - Logical mask indicating which grid points are inside the domain
%       Size: [N x M]
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%	2026/03/05, En-Chi Lee (williameclee@arizona.edu)
%     - Surpressed warnings from POLYSHAPE
%     - Reordered input arguments

function [lon, lat, res, weight, mask] = ...
        polarGridMask(radius, pLonlat, L, options)

    arguments (Input)
        radius (1, 1) {mustBeNumeric, mustBePositive}
        pLonlat (:, 2) {mustBeNumeric}
        L (1, 1) {mustBeNumeric, mustBeInteger} = 18
        options.resFactor (1, 1) {mustBeNumeric, mustBePositive} = 8
    end

    arguments (Output)
        lon (1, :) {mustBeNumeric, mustBeReal}
        lat (:, 1) {mustBeNumeric, mustBeReal}
        res (1, 1) {mustBeNumeric, mustBePositive}
        weight (:, :) {mustBeNumeric, mustBePositive}
        mask (:, :) {mustBeNumericOrLogical}
    end

    % L^2 is the ideal number of points, but we need some buffer since the cap is larger than the region
    numIdealPts = options.resFactor * L ^ 2;
    res = sqrt(360 * radius / numIdealPts);
    res = 360 / ceil(360 / res); % Make sure the resolution can wrap around 360 degrees
    lon = (0:res:360) + res / 2; % Shift by half a resolution to avoid points on the boundary of the grid cells
    lon = lon(1:end - 1); % Remove the last point to avoid duplication at 0 and 360
    lat = (90:-res:(90 - radius)) - res / 2;
    lat = lat(:); % Make sure it's a column vector
    [lonn, latt] = meshgrid(lon, lat);

    weight = cosd(latt) * deg2rad(res) ^ 2;
    % Note: INPOLYGON appears to misclassify points in complex shapes
    warnState = warning("off", "MATLAB:polyshape:repairedBySimplify");
    warnCleanup = onCleanup(@() warning(warnState)); %#ok<NASGU> Ensure warning state is restored
    mask = isinterior(polyshape(pLonlat(:, 1), pLonlat(:, 2)), lonn(:), latt(:));
    mask = reshape(mask, size(lonn));
end
