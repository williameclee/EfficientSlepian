%% DOMAINTOLONLAT - Converts input domain to longitude and latitude coordinates
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function lonlat = domainToLonlat(domain, options)

    arguments (Input)
        domain
        options.AddAnchors (1, 1) {mustBeNumericOrLogical} = false
        options.InputUnit {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
    end

    arguments (Output)
        lonlat (:, 2) {mustBeNumeric}
    end

    if ischar(domain) || isstring(domain)
        lonlat = feval(domain);
    elseif iscell(domain)

        if ~(ischar(domain{1}) || isstring(domain{1}))
            error("If domain is a cell, the first element must be a string or char array representing a function that returns the coordinates of the region boundary.");
        end

        lonlat = feval(domain{:});
    elseif isnumeric(domain)

        if size(domain, 2) ~= 2
            error("If domain is a numeric array, it must have two columns representing the coordinates of the region boundary.");
        end

        lonlat = domain;

        if strcmp(options.InputUnit, 'radians')
            lonlat = rad2deg(lonlat);
        end

    elseif isa(domain, "GeoDomain")
        lonlat = domain.Lonlat("Anchors", options.AddAnchors);
    else
        error("Domain must be a string, cell array, numeric array, or GeoDomain object.");
    end

    % Format the coordinates into a consistent, simple format
    % Is this necessary?
    % [lat, lon] = flatearthpoly(lonlat(:, 2), lonlat(:, 1));
    % lonlat = [lon, lat];

    if any(lonlat(:, 1) < 0)
        warning("Negative longitudes detected. Ensure that the longitude convention is consistent across all inputs and outputs.");
    elseif any(lonlat(:, 1) > 360)
        warning("Longitudes greater than 360 detected. Ensure that the longitude convention is consistent across all inputs and outputs.");
    end

    if strcmp(options.OutputUnit, 'radians')
        lonlat = deg2rad(lonlat);
    end

end
