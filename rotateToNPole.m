%% ROTATETONPOLE - Rotates a domain to the North Pole
%
% Syntax
%   pLonlat = rotateToNPole(domain, pcapLonlat)
%   [pLonlat, pcapLonlat, radius] = rotateToNPole(domain, pcapLonlat)
%
% Input arguments
%   domain - The domain to rotate. Can be a string, cell array, numeric
%       array, or GeoDomain object (see DOMAINTOLONLAT).
%   pcapLonlat - 1x2 vector of the longitude and latitude of the centre of
%       the enclosing cap, in degrees [default: auto-computed]
%
% Output arguments
%   pLonlat - Nx2 array of the longitude and latitude of the rotated
%       domain boundary
%   pcapLonlat - 1x2 vector of the longitude and latitude of the centre
%       of the enclosing cap (in the output unit)
%   radius - scalar radius of the enclosing cap (in the output unit)
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function [pLonlat, pcapLonlat, radius] = rotateToNPole(domain, pcapLonlat, options)

    arguments (Input)
        domain
        pcapLonlat (1, 2) {mustBeNumeric, mustBeReal} = [-999, -999]
        options.InputUnit {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
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

    lonlat = domainToLonlat(domain, "AddAnchors", true, "InputUnit", options.InputUnit, "OutputUnit", "radians");
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
        % pcapLonlat is already in radians
    end

    if nargout > 0
        return
    end

    %% Visualisation
    figure
    plot(lonlat(:, 2), lonlat(:, 1), 'k')

    axis equal tight
end
