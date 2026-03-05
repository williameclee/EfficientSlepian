%% ENCLOSINGCAP - Finds the centre and radius of the smallest enclosing spherical cap
% Syntax
%   [lonlat, radius] = enclosingCap(domain)
%
% Input arguments
%     - A string or char: name of a function returning boundary coordinates
%       For example, "antarctica" (from slepian_delta) or "npacific" (from
%       ULMO)
%     - A cell array: {funcName, args...} passed to feval
%       Like the above, funcName is the name of a function that returns
%       boundary coordinates, and args are additional arguments to that
%       function.
%     - A numeric Nx2 array: longitude-latitude boundary coordinates
%       The unit is specified by the InputUnit option.
%     - A GeoDomain object (from the ULMO package)
%       See DOMAONTOLONLAT for details.
%   InputUnit (name-value) - Unit of the input coordinates ("degrees" or
%       "radians")
%       This option only matters if the input domain is a numeric array.
%       The default unit is "degrees".
%   OutputUnit (name-value) - Unit of the output coordinates ("degrees" or
%       "radians")
%       The default unit is "degrees".
%
% Output arguments
%   lonlat - Longitude and latitude of the centre of the enclosing cap, in
%       the output unit.
%   radius - Radius of the enclosing cap, in the output unit.
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function [pcapLonlat, radius] = enclosingCap(domain, options)

    arguments (Input)
        domain
        options.InputUnit {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
    end

    arguments (Output)
        pcapLonlat (1, 2) {mustBeNumeric}
        radius (1, 1) {mustBeNumeric, mustBePositive}
    end

    lonlat = domainToLonlat(domain, "AddAnchors", true, "InputUnit", options.InputUnit, "OutputUnit", "degrees");

    % Find the longest distance between 2 points on the boundary
    maxDists = zeros(length(lonlat), 1);
    maxDistLonlats = zeros(length(lonlat), 4);

    for i = 1:length(lonlat)
        dists = distance(lonlat(i, 2), lonlat(i, 1), lonlat(:, 2), lonlat(:, 1));
        [maxDists(i), idx] = max(dists);
        maxDistLonlats(i, :) = [lonlat(i, :), lonlat(idx, :)];
    end

    [maxDist, idx] = max(maxDists);
    maxDistLonlat = maxDistLonlats(idx, :);

    % Find the centre of the cap as the midpoint of the longest distance
    p1Xyz = [0, 0, 0];
    p2Xyz = [0, 0, 0];
    [p1Xyz(1), p1Xyz(2), p1Xyz(3)] = sph2cart( ...
        deg2rad(maxDistLonlat(1)), deg2rad(maxDistLonlat(2)), 1);
    [p2Xyz(1), p2Xyz(2), p2Xyz(3)] = sph2cart( ...
        deg2rad(maxDistLonlat(3)), deg2rad(maxDistLonlat(4)), 1);
    pmXyz = (p1Xyz + p2Xyz) / 2;
    [pcapLonr, pcapLatr, ~] = cart2sph(pmXyz(1), pmXyz(2), pmXyz(3));
    pcapLonlat = [wrapTo360(rad2deg(pcapLonr)), rad2deg(pcapLatr)];

    radius = maxDist / 2;

    if strcmp(options.OutputUnit, 'radians')
        pcapLonlat = deg2rad(pcapLonlat);
        radius = deg2rad(radius);
    end

    if nargout > 0
        return
    end

    %% Visualisation
    figure
    plot(lonlat(:, 1), lonlat(:, 2), 'k');
    hold on
    scatter(maxDistLonlat([1, 3]), maxDistLonlat([2, 4]), 'b')
    plot([maxDistLonlat(1), pcapLonlat(1), maxDistLonlat(3)], ...
        [maxDistLonlat(2), pcapLonlat(2), maxDistLonlat(4)], ...
    'b')
    scatter(pcapLonlat(1), pcapLonlat(2), 'r')
    hold off

    title('Enclosing Cap')
    xlabel('Longitude')
    ylabel('Latitude')
    axis equal tight
end
