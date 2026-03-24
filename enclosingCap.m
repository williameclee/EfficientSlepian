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
%   maxIters (name-value) - Maximum number of iterations for the algorithm
%       The default number is 1000.
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
%       in the output unit
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

    if nargout > 0

        if strcmp(options.OutputUnit, 'radians')
            clonlat = deg2rad(clonlat);
            radius = deg2rad(radius);
        end

        return
    end

    %% Visualisation
    makeDemoPlot(lonlat, clonlat, radius)
end

%% Subfunctions
function makeDemoPlot(lonlat, clonlat, radius)
    xyz = [cosd(lonlat(:, 2)) .* cosd(lonlat(:, 1)), ...
               cosd(lonlat(:, 2)) .* sind(lonlat(:, 1)), ...
               sind(lonlat(:, 2))];

    % The boundary of the polar cap in Cartesian coordinates
    nPts = 100;
    pcapXyz = [zeros(nPts, 2), ones(nPts, 1)] * cosd(radius);
    pcapXyz(:, 1) = sind(radius) * cosd(linspace(0, 360, nPts));
    pcapXyz(:, 2) = sind(radius) * sind(linspace(0, 360, nPts));

    % Rotate the boundary points to be centered at capLonlat
    theta = 90 - clonlat(2); % polar angle
    phi = clonlat(1); % azimuthal angle
    Rcol = [cosd(theta), 0, sind(theta); 0, 1, 0; -sind(theta), 0, cosd(theta)];
    Rlon = [cosd(phi), sind(phi), 0; sind(phi), -cosd(phi), 0; 0, 0, 1];
    R = Rlon * Rcol;
    capXyz = (R * pcapXyz')';
    capLonlat = ...
        [atan2d(capXyz(:, 2), capXyz(:, 1)), ...
         asind(capXyz(:, 3))];
    capLonlat(:, 1) = mod(capLonlat(:, 1), 360); % Ensure longitudes are in [0, 360]

    % Visualize the points and the enclosing cap
    figure
    subplot(1, 2, 1)

    hold on
    plot(lonlat(:, 1), lonlat(:, 2), 'b')
    plot(capLonlat(:, 1), capLonlat(:, 2), 'r')
    scatter(clonlat(1), clonlat(2), 'ro', 'filled')
    hold off
    axis equal tight
    xlabel('Longitude [°]')
    ylabel('Latitude [°]')

    subplot(1, 2, 2)
    hold on
    plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'b')
    plot3(capXyz(:, 1), capXyz(:, 2), capXyz(:, 3), 'r')
    scatter3(cosd(clonlat(2)) * cosd(clonlat(1)), cosd(clonlat(2)) * sind(clonlat(1)), sind(clonlat(2)), 'ro', 'filled')
    hold off
    axis equal tight

    % Change the view angle for better visualization
    view(90 + clonlat(1), clonlat(2))

    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end
