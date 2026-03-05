%% GLMALPHA_EFF - Computes the Slepian basis using the efficient way
%
% Syntax
%   [G, V, N] = glmalpha_eff(domain, L)
%   [G, V, N] = glmalpha_eff(domain, L, "Name", value)
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
%       See DOMAINTOLONLAT for details.
%   L (optional) - Bandwidth (maximum angular degree)
%       The default degree is 18.
%   pcapConcThreshold (name-value) - Minimum energy concentration value for
%       a polar-cap Slepian function to not be discarded.
%       The default value is 0.3.
%
% Output arguments
%   G - Projection matrix from the Slepian basis to the spherical harmonics
%       Size: [(L+1)^2 x numFuns]
%   V - Concentration eigenvalues in descending order
%       Size: [numFuns x 1]
%   N - Shannon number
%       Estimated number of well-concentrated functions, proportional to the
%       area of the domain and the squared bandwidth.
%
% See also
%   GLMALPHA, GRUNBAUM
%
% Author
%	2026/03/05, En-Chi Lee (williameclee@arizona.edu)

function [G, V, N] = glmalpha_eff(domain, L, options)

    arguments (Input)
        domain
        L (1, 1) {mustBeInteger, mustBePositive} = 18
        options.pcapConcThreshold (1, 1) {mustBePositive} = 0.3
    end

    arguments (Output)
        G (:, :) {mustBeReal}
        V (:, 1) {mustBeNonnegative}
        N (1, 1) {mustBePositive}
    end

    % Step 1: Find the enclosing polar cap
    [pcapLonlatd, radiusd] = enclosingCap(domain, "OutputUnit", "degrees");

    % Step 2: Rotate the domain to the North Pole
    pLonlatd = rotateToNPole(domain, pcapLonlatd, "OutputUnit", "degrees");

    % Step 3: compute the Slepian functions for the polar cap
    pcapConcThreshold = min(options.pcapConcThreshold, 1); % Threshold for well-concentrated Slepian functions

    % Preparation for Step 3c: Make the colatitude and longitude grid for evaluating the Slepian functions spatially
    [pgridLond, pgridLatd, ~, pgridWeight, pgridMask] = ...
        polarGridMask(radiusd, pLonlatd, L, resFactor = 16);

    pcapConcs = []; % The eigenvalues
    pcapGs = {}; % The Slepian coefficients
    pcapMs = []; % The order each function corresponds to
    pcapSlepMesh = []; % The Slepian functions evaluated on the grid

    for m = -L:L
        % Step 3a: Find the spherical harmonic coefficients of Slepian functions of order m
        [~, ~, ~, pcapG_m, ~, pcapConc_m] = grunbaum_new(radiusd, L, m, 0);

        % Step 3b: Discard poorly concentrated Slepian functions
        numConc = sum(pcapConc_m > pcapConcThreshold);

        if numConc == 0
            continue
        end

        pcapG_m = pcapG_m(:, 1:numConc);
        pcapConc_m = pcapConc_m(1:numConc);
        % Don't care about this yet
        % G = [G; G_m(:)];
        pcapConcs = [pcapConcs; pcapConc_m(:)];
        pcapMs = [pcapMs; repmat(m, numConc, 1)];

        for i = 1:numConc
            pcapGs = [pcapGs; {pcapG_m(:, i)}];
        end

        % Step 3c: Evaluate the Slepian functions spatially and store in slep
        Ylm_m = zeros(length(pgridLatd), length(pgridLond), L - abs(m) + 1);

        for l = abs(m):L
            Ylm_m(:, :, l - abs(m) + 1) = ylm(l, m, deg2rad(90 - pgridLatd), deg2rad(pgridLond));
        end

        for i = 1:numConc
            pcapSlepMesh_i = sum(Ylm_m .* reshape(pcapG_m(:, i), 1, 1, []), 3);
            pcapSlepMesh = cat(3, pcapSlepMesh, pcapSlepMesh_i);
        end

    end

    % Sort the Slepian functions by concentration
    [pcapConcs, pcapConcSortId] = sort(pcapConcs, 'descend');
    pcapSlepMesh = pcapSlepMesh(:, :, pcapConcSortId);
    pcapGs = pcapGs(pcapConcSortId);
    pcapMs = pcapMs(pcapConcSortId);
    numFuns = length(pcapConcs);

    % Step 4: Compute localisation matrix P
    locMat = nan(numFuns, numFuns);

    for i = 1:numFuns
        locMat(i, i) = sum(pcapSlepMesh(:, :, i) .^ 2 .* pgridMask .* pgridWeight, "all");

        for j = i + 1:numFuns
            locMat(i, j) = sum( ...
                pcapSlepMesh(:, :, i) .* pcapSlepMesh(:, :, j) .* pgridMask .* pgridWeight, "all");
            locMat(j, i) = locMat(i, j);
        end

    end

    % Step 5: Eigen-decomposition of localisation matrix
    [locMatEigvecs, locMatEigvals] = eig(locMat);
    [locMatEigvals, eigSortId] = sort(diag(locMatEigvals), 'descend');
    locMatEigvecs = locMatEigvecs(:, eigSortId);

    % Step 6: Get the Slepian functions for the rotated domain
    pG = zeros((L + 1) ^ 2, numFuns); % Same as GLMALPHA

    for i = 1:numFuns
        m = pcapMs(i);
        pcapG_m = pcapGs{i};
        pG((abs(m):L) .* ((abs(m):L) + 1) + m + 1, i) = pcapG_m;
    end

    pG = pG * locMatEigvecs;

    % Rotate the Slepian functions back to the original domain
    G = rotateG(L, pG, pcapLonlatd);

    V = locMatEigvals(:); % Concentration ratios (eigenvalues)
    N = (L + 1) ^ 2 * spharea(pLonlatd);
end

%% Subfunctions
function G = rotateG(L, pG, pcapLonlatd)
    [degrees, orders, ~, lmcosi, ~, mzo, ~, ~, rinm, ronm, ~] = addmon(L);

    numFuns = size(pG, 2);
    CC = cell([1, numFuns]);

    parfor j = 1:size(pG, 2)
        cosi = lmcosi(:, 3:4); % Blank coefficient template
        cosi(ronm) = pG(:, j); % Insert North Pole coefficients
        CC{j} = cosi;
    end

    CC_coeff = cell(1, numFuns);

    parfor i = 1:numFuns
        CC_coeff{i} = kindeks( ...
            plm2rot([orders, degrees, CC{i}], 180, - (90 - pcapLonlatd(2)), -pcapLonlatd(1)), 3:4);
    end

    G = nan((L + 1) ^ 2, numFuns);

    parfor j = 1:size(pG, 2)
        cosi = CC_coeff{j};
        % Remove the m=0 sine coefficients (which are always zero)
        cosinozero = cosi(mzo);
        % Reorder into standard lmcosi format
        G(:, j) = cosinozero(rinm);
    end

end
