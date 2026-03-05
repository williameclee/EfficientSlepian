%% ROTATETONPOLE - Rotates a domain to the North Pole
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
    elseif strcmp(options.OutputUnit, "radians")
        pLonlat = deg2rad(pLonlat);
        pcapLonlat = deg2rad(pcapLonlat);
    end

    if nargout > 0
        return
    end

    %% Visualisation
    figure
    plot(lonlat(:, 2), lonlat(:, 1), 'k')

    axis equal tight
end
