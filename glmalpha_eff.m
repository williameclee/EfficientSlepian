%% GLMALPHA_EFF - Computes the Slepian basis using the efficient way
% We follow the Formulation of Bates, Alice P et al. (2017) to compute the
% Slepian basis for an arbitrary domain. The 'Steps' 1 to 6 mentioned below
% map one-to-one to the steps in Section III-E of the paper (p. 4385).
% In this efficient formulation, there are essentially two stages to
% compute the Slepian basis:
%    1. Compute the Slepian functions for the polar cap that encloses the
%       domain, which forms a much smaller basis.
%    2. Compute, again, the Slepian functions for the polar cap Slepian
%       basis over the rotated domain.
% The 'Slepian functions of the Slepian functions' are then projected back
% to the spherical harmonics and rotated back to the original domain to get
% the final projection matrix G.
% Because the first stage can exploit the axisymmetry of the polar cap, and
% the second stage is only computed over a much smaller number of basis
% functions, this efficient formulation can be much faster than the direct
% formulation (GLMALPHA) for large bandwidths and small domains.
%
% Syntax
%   [G, V, N] = glmalpha_eff(domain, L)
%   [G, V, N] = glmalpha_eff(domain, L, truncation, rotb)
%   [G, V, N] = glmalpha_eff(__, "Name", value)
%
% Input arguments
%   domain - Domain to convert
%     - A string or char: name of a function returning boundary coordinates
%       For example, "antarctica" (from slepian_delta) or "npacific" (from
%       ULMO)
%     - A cell array: {funcName, args...} passed to FEVAL
%       Like the above, funcName is the name of a function that returns
%       boundary coordinates, and args are additional arguments to that
%       function.
%     - A numeric Nx2 array: longitude-latitude boundary coordinates
%       The unit is specified by the InputUnit option.
%     - A GeoDomain object (from the ULMO package)
%       See DOMAINTOLONLAT for details.
%   L (optional) - Bandwidth (maximum angular degree)
%       The default degree is 18.
%   truncation (optional) - Number of final Slepian functions to return, in
%       descending order of concentration
%       If not specified, all computed functions will be returned; if set
%       to "N" (case-insensitive), the number of functions will be
%       determined by the Shannon number (N).
%       The default truncation level is not specified (i.e. no truncation).
%   rotb (optional) - Whether to rotate the Slepian functions back to the
%       original domain
%       If false, the returned projection matrix will be for the domain
%       rotated to the north pole.
%       The default value is TRUE.
%   pcapConcThreshold (name-value) - Minimum energy concentration value for
%       a polar-cap Slepian function to not be discarded.
%       The default value is 0.3.
%   resFactor (name-value) - Resolution factor for the polar grid used to
%       compute the second localisation matrix (of the polar cap Slepian
%       basis)
%       The integral for the localisation matrix will be evaluated over a
%       polar grid with roughly L^2 * resFactor points. A higher resolution
%       factor should give a more accurate localisation matrix, but will
%       take more time to compute.
%       The default value is 8.
%
% Output arguments
%   G - Projection matrix from the Slepian basis to the spherical harmonics
%       Size: [(L+1)^2 x numFuns], where numFuns is the number of Slepian
%       functions returned (after truncation, if applicable)
%   V - Concentration eigenvalues in descending order
%       The first row contains the eigenvalues of the polar cap Slepian
%       functions, and the second row contains the eigenvalues of the
%       Slepian functions for the rotated domain relative to the polar cap
%       Slepian basis. That is, they are not comparable to the eigenvalues
%       from GLMALPHA.
%       Size: [2 x numFuns]
%   N - Shannon number
%       Estimated number of well-concentrated functions, proportional to the
%       area of the domain and the squared bandwidth.
%
% See also
%   GLMALPHA, GRUNBAUM
%
% Author
%	2026/03/05, En-Chi Lee (williameclee@arizona.edu)
% Last modified
%	2026/03/06, En-Chi Lee (williameclee@arizona.edu)
%     - Changed eigenvalues (V) output format
%     - Added truncation and rotb arguments

function [G, V, N] = glmalpha_eff(domain, L, truncation, rotb, options)

    arguments (Input)
        domain
        L (1, 1) {mustBeInteger, mustBePositive} = 18
        truncation (1, 1) = NaN
        rotb (1, 1) {mustBeNumericOrLogical} = true
        options.pcapConcThreshold (1, 1) ...
            {mustBeInRange(options.pcapConcThreshold, 0, 1, "exclude-lower")} = 0.3
        options.resFactor (1, 1) {mustBePositive} = 8
    end

    arguments (Output)
        G (:, :) {mustBeReal}
        V (2, :) {mustBeReal, mustBeInRange(V, 0, 1)}
        N (1, 1) {mustBePositive}
    end

    if isnumeric(truncation)

        if ~isnan(truncation) && (truncation <= 0)
            error("Truncation must be positive, but got %s.", num2str(truncation));
        end

    elseif (isstring(truncation) || ischar(truncation))

        if ~strcmpi(truncation, "N")
            error("Truncation must be either a positive value or the string 'N', but got '%s'.", truncation);
        end

    else
        error("Truncation must be either a positive value or the string 'N', but got class %s.", ...
            upper(class(truncation)));
    end

    pcapConcThreshold = options.pcapConcThreshold;
    resFactor = options.resFactor;

    %% Main computation
    % Step 1: Find the enclosing polar cap
    [pcapLonlatd, radiusd] = enclosingCap(domain, "OutputUnit", "degrees");

    % Step 2: Rotate the domain to the North Pole
    pLonlatd = rotateToNPole(domain, pcapLonlatd, "OutputUnit", "degrees");

    % Step 3: compute the Slepian functions for the polar cap
    % Preparation for Step 3c: Make the colatitude and longitude grid for
    % evaluating the Slepian functions spatially
    [pgridLond, pgridLatd, ~, pgridWeight, pgridMask] = ...
        polarGridMask(radiusd, pLonlatd, L, resFactor = resFactor);
    pgridWeight = pgridWeight .* pgridMask; % Mask the weights

    pcapConcs = []; % The eigenvalues
    pcapGs = {}; % The Slepian coefficients
    pcapMs = []; % The order each function corresponds to
    pcapSlepMesh = []; % The Slepian functions evaluated on the grid

    for m = -L:L
        % Step 3a: Find the SH coefficients of polar cap Slepian functions
        % of order m
        [~, ~, ~, pcapG_m, ~, pcapConc_m] = grunbaum_new(radiusd, L, m, 0);

        % Step 3b: Discard poorly concentrated Slepian functions
        numConc = sum(pcapConc_m > pcapConcThreshold);

        if numConc == 0
            continue
        end

        pcapG_m = pcapG_m(:, 1:numConc);
        pcapConc_m = pcapConc_m(1:numConc);
        pcapConcs = [pcapConcs; pcapConc_m(:)];
        pcapMs = [pcapMs; repmat(m, numConc, 1)];

        for i = 1:numConc
            pcapGs = [pcapGs; {pcapG_m(:, i)}];
        end

        % Step 3c: Evaluate the Slepian functions spatially and store in slep
        Ylm_m = zeros(length(pgridLatd), length(pgridLond), L - abs(m) + 1);

        for l = abs(m):L
            Ylm_m(:, :, l - abs(m) + 1) = ...
                ylm(l, m, deg2rad(90 - pgridLatd), deg2rad(pgridLond));
        end

        for i = 1:numConc
            pcapSlepMesh_i = sum(Ylm_m .* reshape(pcapG_m(:, i), 1, 1, []), 3);
            pcapSlepMesh = cat(3, pcapSlepMesh, pcapSlepMesh_i);
        end

    end

    % Sort the Slepian functions by concentration
    [pcapConcs, pcapConcSortId] = sort(pcapConcs, "descend");
    pcapSlepMesh = pcapSlepMesh(:, :, pcapConcSortId);
    pcapGs = pcapGs(pcapConcSortId);
    pcapMs = pcapMs(pcapConcSortId);
    numFuns = length(pcapConcs);

    % Projection matrix from the polar cap Slepian basis to SH coefficients
    pcapG = zeros((L + 1) ^ 2, numFuns);

    for i = 1:numFuns
        m = pcapMs(i);
        pcapG_m = pcapGs{i};
        pcapG((abs(m):L) .* ((abs(m):L) + 1) + m + 1, i) = pcapG_m;
    end

    % Step 4: Compute localisation matrix for the polar cap Slepian basis
    % over the rotated domain
    locMat = nan(numFuns, numFuns);

    for i = 1:numFuns
        locMat(i, i) = sum(pcapSlepMesh(:, :, i) .^ 2 .* pgridWeight, "all");

        for j = i + 1:numFuns
            locMat(i, j) = sum( ...
                pcapSlepMesh(:, :, i) .* pcapSlepMesh(:, :, j) .* pgridWeight, "all");
            locMat(j, i) = locMat(i, j);
        end

    end

    % Step 5: Eigen-decomposition of localisation matrix
    % Get the Slepian functions for the polar cap Slepian functions
    [pSlepG, pSlepConcs] = eig(locMat);
    [pSlepConcs, pConcSortId] = sort(diag(pSlepConcs), "descend");
    pSlepG = pSlepG(:, pConcSortId);

    % Step 6: Get the Slepian functions for the rotated domain
    % Project back to the spherical harmonics
    pG = pcapG * pSlepG;

    % Rotate the Slepian functions back to the original domain
    if rotb
        G = rotateG(L, pG, pcapLonlatd);
    else
        G = pG;
    end

    V = [pcapConcs(:), pSlepConcs(:)].'; % Concentrations (eigenvalues)
    N = (L + 1) ^ 2 * spharea(pLonlatd);

    % Truncate the basis if asked
    if (isstring(truncation) || ischar(truncation)) && strcmpi(truncation, "N")
        truncation = round(N);
    end

    if ~isnan(truncation)
        truncation = round(truncation); % In case it's a non-integer numeric value

        if truncation <= numFuns
            G = G(:, 1:truncation);
            V = V(:, 1:truncation);
        else
            warning( ...
                "Truncation level (%d) is larger than the number of computed functions (%d). No truncation applied.", truncation, numFuns);
        end

    end

end

%% Subfunctions
function G = rotateG(L, pG, pcapLonlatd)
    % Rotates the projection matrix (pG) from the north pole back to the
    % original domain (G)

    arguments (Input)
        L (1, 1) {mustBeInteger, mustBePositive}
        pG (:, :) {mustBeReal}
        pcapLonlatd (1, 2) {mustBeReal}
    end

    arguments (Output)
        G (:, :) {mustBeReal}
    end

    [degrees, orders, ~, lmcosi, ~, mzo, ~, ~, rinm, ronm] = addmon(L);

    numFuns = size(pG, 2);
    CC = cell([1, numFuns]);

    parfor j = 1:size(pG, 2)
        cosi = lmcosi(:, 3:4); % Blank coefficient template
        cosi(ronm) = pG(:, j); % Insert North Pole coefficients
        CC{j} = cosi;
    end

    CC_coeff = cell(1, numFuns);

    rotLatd = (90 - pcapLonlatd(2));
    rotLond = pcapLonlatd(1);

    parfor i = 1:numFuns
        CC_coeff{i} = kindeks( ...
            plm2rot([orders, degrees, CC{i}], 180, -rotLatd, -rotLond), 3:4);
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
