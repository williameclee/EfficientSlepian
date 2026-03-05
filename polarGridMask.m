%% POLARGRIDMASK - Creates a grid in the polar cap and a corresponding domain mask
%
% Syntax
%   [lon, lat, res, weight, mask] = polarGridMask(radius, L, pLonlat)
%   [lon, lat, res, weight, mask] = polarGridMask(radius, L, pLonlat, resFactor=N)
%
% Input arguments
%   radius - scalar radius of the polar cap, in degrees
%   L - bandwidth (maximum angular degree), used to determine grid resolution
%   pLonlat - Nx2 array of the rotated domain boundary [longitude, latitude]
%       in degrees
%   resFactor - resolution scaling factor; higher values produce finer
%       grids [default: 8]
%
% Output arguments
%   lon - 1xM vector of grid longitudes, in degrees
%   lat - Nx1 vector of grid latitudes, in degrees
%   res - scalar grid resolution, in degrees
%   weight - NxM matrix of integration weights (solid-angle element)
%   mask - NxM logical matrix, true inside the rotated domain
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function [lon, lat, res, weight, mask] = ...
        polarGridMask(radius, L, pLonlat, options)

    arguments (Input)
        radius (1, 1) {mustBeNumeric, mustBePositive}
        L (1, 1) {mustBeNumeric, mustBeInteger}
        pLonlat (:, 2) {mustBeNumeric}
        options.resFactor (1, 1) {mustBeNumeric, mustBePositive} = 8
    end

    arguments (Output)
        lon (1, :) {mustBeNumeric, mustBeReal}
        lat (:, 1) {mustBeNumeric, mustBeReal}
        res (1, 1) {mustBeNumeric, mustBePositive}
        weight (:, :) {mustBeNumeric, mustBePositive}
        mask (:, :) {mustBeNumericOrLogical}
    end

    numIdealPts = options.resFactor * L ^ 2; % L^2 is the ideal number of points, but we need some buffer since the cap is larger than the region
    res = sqrt(360 * radius / numIdealPts);
    res = 360 / ceil(360 / res); % Make sure the resolution can wrap around 360 degrees
    lon = (0:res:360) + res / 2; % Shift by half a resolution to avoid points on the boundary of the grid cells
    lon = lon(1:end - 1); % Remove the last point to avoid duplication at 0 and 360
    lat = (90:-res:(90 - radius)) - res / 2;
    lat = lat(:); % Make sure it's a column vector
    [lonn, latt] = meshgrid(lon, lat);

    weight = cosd(latt) * deg2rad(res) ^ 2;
    % Note: INPOLYGON appears to misclassify points in complex shapes
    mask = isinterior(polyshape(pLonlat(:, 1), pLonlat(:, 2)), lonn(:), latt(:));
    mask = reshape(mask, size(lonn));
end
