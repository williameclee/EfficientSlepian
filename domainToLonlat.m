%% DOMAINTOLONLAT - Converts input domain to longitude-latitude coordinates
%
% Syntax
%   lonlat = domainToLonlat(domain)
%
% Input arguments
%   domain - Domain to convert
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
%   AddAnchors (name-value) - Whether to add anchor points to the boundary
%       Adding anchor points along long straight boundaries can help make
%       sure the shape stays intact when rotated.
%       The default behaviour is FALSE.
%   InputUnit (name-value) - Unit of the input coordinates ("degrees" or
%       "radians")
%       This option only matters if the input domain is a numeric array.
%       The default unit is "degrees".
%   OutputUnit (name-value) - Unit of the output coordinates ("degrees" or
%       "radians")
%       The default unit is "degrees".
%
% Output arguments
%   lonlat - Nx2 array of boundary coordinates [longitude, latitude]
%
% See also
%   GEODOMAIN (from ULMO)
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)

function lonlat = domainToLonlat(domain, options)

    arguments (Input)
        domain
        options.AddAnchors (1, 1) {mustBeNumericOrLogical} = false
        options.InputUnit ...
            {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit ...
            {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
    end

    arguments (Output)
        lonlat (:, 2) {mustBeNumeric}
    end

    if ischar(domain) || isstring(domain)
        lonlat = feval(domain);
    elseif iscell(domain)

        if ~(ischar(domain{1}) || isstring(domain{1}))
            error("Slepian:Domain:InvalidInputType", ...
                "If domain is a cell, the first element must be the name of a function, but got %s type.", ...
                upper(class(domain{1})));
        end

        lonlat = feval(domain{:});
    elseif isnumeric(domain)

        if ndims(domain) ~= 2 || size(domain, 2) ~= 2
            error("Slepian:Domain:InvalidInputType", ...
                "If domain is a numeric array, it must have two columns representing the coordinates of the region boundary, but got an array of size %s.", ...
                mat2str(size(domain)));
        end

        lonlat = domain;

        if strcmp(options.InputUnit, 'radians')
            lonlat = rad2deg(lonlat);
        end

    elseif isa(domain, "GeoDomain")
        lonlat = domain.Lonlat("Anchors", options.AddAnchors);
    else
        error("Slepian:Domain:InvalidInputType", ...
            "Domain must be a string, cell array, numeric array, or GeoDomain object, but got %s type.", ...
            upper(class(domain)));
    end

    if strcmp(options.OutputUnit, 'radians')
        lonlat = deg2rad(lonlat);
    end

end
