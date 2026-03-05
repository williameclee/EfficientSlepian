%% DOMAINTOLONLAT - Converts input domain to longitude and latitude coordinates
%
% Syntax
%   lonlat = domainToLonlat(domain)
%
% Input arguments
%   domain - The domain to convert. Can be:
%       - A string or char: name of a function returning boundary coordinates
%       - A cell array: {funcName, args...} passed to feval
%       - A numeric Nx2 array: longitude-latitude boundary coordinates
%       - A GeoDomain object
%   AddAnchors - Whether to add anchor points to the boundary [default: false]
%   InputUnit - Unit of the input coordinates ('degrees' or 'radians')
%       [default: 'degrees']
%   OutputUnit - Unit of the output coordinates ('degrees' or 'radians')
%       [default: 'degrees']
%
% Output arguments
%   lonlat - Nx2 array of boundary coordinates [longitude, latitude]
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

    if any(lonlat(:, 1) < 0)
        warning("Negative longitudes detected. Ensure that the longitude convention is consistent across all inputs and outputs.");
    end

    if any(lonlat(:, 1) > 360)
        warning("Longitudes greater than 360 detected. Ensure that the longitude convention is consistent across all inputs and outputs.");
    end

    if strcmp(options.OutputUnit, 'radians')
        lonlat = deg2rad(lonlat);
    end

end
