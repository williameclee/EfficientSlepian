%% ROTATETONPOLE - Rotates a domain to the North Pole
%
% Syntax
%   pLonlat = rotateToNPole(domain, pcapLonlat)
%   [pLonlat, pcapLonlat, radius] = rotateToNPole(domain, pcapLonlat)
%
% Input arguments
%   domain - Domain to rotate
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
%   pcapLonlat (optional) - 1x2 vector of the longitude and latitude of the
%       centre of the enclosing cap, in degrees
%       If not provided, the function will compute the enclosing cap and
%       use its centre.
%   InputUnit (name-value) - Unit of the input coordinates ("degrees" or
%       "radians")
%       This option only matters if the input domain is a numeric array.
%       The default unit is "degrees".
%   OutputUnit (name-value) - Unit of the output coordinates ("degrees" or
%       "radians")
%       The default unit is "degrees".
%
% Output arguments
%   pLonlat - [longitude, latitude] coordinates of the rotated domain
%       boundary
%       Size: [N x 2]
%   pcapLonlat - [longitude, latitude] coordinates of the centre of the
%       enclosing cap
%       Size: [1 x 2]
%   radius - Radius of the enclosing cap (in the output unit)
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function [pLonlat, pcapLonlat, radius] = rotateToNPole(domain, pcapLonlat, options)

    arguments (Input)
        domain
        pcapLonlat (1, 2) {mustBeNumeric, mustBeReal} = [-999, -999]
        options.InputUnit ...
            {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit ...
            {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
    end

    arguments (Output)
        pLonlat (:, 2) {mustBeNumeric}
        pcapLonlat (1, 2) {mustBeNumeric, mustBeReal}
        radius (1, 1) {mustBeNumeric}
    end

    if isequal(pcapLonlat, [-999, -999])
        [pcapLonlat, radius] = enclosingCap(domain, ...
            "InputUnit", options.InputUnit, "OutputUnit", "radians");
    elseif strcmp(options.InputUnit, 'degrees')
        pcapLonlat = deg2rad(pcapLonlat);
        radius = NaN; % Not calculated
    end

    pcapCollon = [pi / 2 - pcapLonlat(2), pcapLonlat(1)];

    lonlat = domainToLonlat(domain, "AddAnchors", true, ...
        "InputUnit", options.InputUnit, "OutputUnit", "radians");
    collon = [pi / 2 - lonlat(:, 2), lonlat(:, 1)];

    %% Rotation to 0° longitude and 90° latitude (North Pole)
    [pCol, pLon] = rottp( ...
        collon(:, 1), collon(:, 2), pcapCollon(2), pcapCollon(1), 0);
    pLatd = 90 - rad2deg(pCol);
    pLond = rad2deg(pLon);
    [pLatd, pLond] = flatearthpoly(pLatd, pLond, 180);
    pLonlat = [pLond, pLatd];

    if strcmp(options.OutputUnit, "degrees")
        radius = rad2deg(radius);
        pcapLonlat = rad2deg(pcapLonlat);
    elseif strcmp(options.OutputUnit, "radians")
        pLonlat = deg2rad(pLonlat);
    end

    if nargout > 0
        return
    end

    %% Visualisation
    figure
    plot(lonlat(:, 2), lonlat(:, 1), 'k')

    axis equal tight
end
