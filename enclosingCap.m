%% ENCLOSINGCAP - Finds smallest enclosing spherical cap of shape on sphere
% The algorithm starts from the centroid of the input points and
% iteratively moves towards the furthest point until convergence. As an
% extra final step, the radius is calculated to ensure all points are
% enclosed.
%
% Syntax
%   [clonlat, radius] = enclosingCap(domain)
%   [clonlat, radius] = enclosingCap(domain, "Name", value)
%
% Input arguments
%   domain - Domain to be enclosed by the cap
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
%       See DOMAINTOLONLAT for details.
%   InputUnit (name-value) - Unit of the input coordinates ("degrees" or
%       "radians")
%       This option only matters if the input domain is a numeric array.
%       The default unit is "degrees".
%   OutputUnit (name-value) - Unit of the output coordinates ("degrees" or
%       "radians")
%       The default unit is "degrees".
%
% Output arguments
%   clonlat - Longitude and latitude of the centre of the enclosing cap,
%       in the output unit.
%       Size: [1 x 2]
%   radius - Radius of the enclosing cap, in the output unit
%
% Author
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%	2026/03/23, En-Chi Lee (williameclee@arizona.edu)
%     - Switched to iterative algorithm with the correct objective

function [clonlat, radius] = enclosingCap(domain, options)

    arguments (Input)
        domain
        options.maxIters (1, 1) {mustBePositive, mustBeInteger} = 1000
        options.InputUnit ...
            {mustBeMember(options.InputUnit, {'degrees', 'radians'})} = 'degrees'
        options.OutputUnit ...
            {mustBeMember(options.OutputUnit, {'degrees', 'radians'})} = 'degrees'
    end

    arguments (Output)
        clonlat (1, 2) {mustBeNumeric, mustBeFinite}
        radius (1, 1) {mustBeNumeric, mustBePositive}
    end

    lonlat = domainToLonlat(domain, "AddAnchors", true, ...
        "InputUnit", options.InputUnit, "OutputUnit", "degrees");

    lonlat = lonlat(~any(isnan(lonlat), 2), :); % Remove any rows with NaN values

    xyz = ...
        [cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
         cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
         sind(lonlat(:, 2))];

    % Start from the centroid
    cxyz = mean(xyz, 1);
    cxyz = cxyz / norm(cxyz);

    for iIter = 1:options.maxIters
        % Find the furthest point
        dists = xyz * cxyz';
        [~, jMin] = min(dists);
        farthestXyz = xyz(jMin, :);

        % Move the centre towards the furthest point
        stepSize = 1 / (iIter + 1);
        cxyz = cxyz + stepSize * (farthestXyz - cxyz);
        cxyz = cxyz / norm(cxyz);
    end

    % Project back to lonlat and calculate the radius
    clonlat = [wrapTo360(atan2d(cxyz(2), cxyz(1))), asind(cxyz(3))];
    radius = acosd(min(min(xyz * cxyz'), 1));

    if strcmp(options.OutputUnit, 'radians')
        clonlat = deg2rad(clonlat);
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
    plot([maxDistLonlat(1), clonlat(1), maxDistLonlat(3)], ...
        [maxDistLonlat(2), clonlat(2), maxDistLonlat(4)], ...
    'b')
    scatter(clonlat(1), clonlat(2), 'r')
    hold off

    title('Enclosing Cap')
    xlabel('Longitude')
    ylabel('Latitude')
    axis equal tight
end
