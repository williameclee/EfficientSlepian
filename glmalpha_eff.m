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
%   [G, V, N, K] = glmalpha_eff(domain, L)
%   [G, V, N, K] = glmalpha_eff(domain, L, truncation, rotb)
%   [G, V, N, K] = glmalpha_eff(__, "Name", value)
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
%   IntegrationMethod (name-value) - Method to compute the localisation
%       matrix for the polar cap Slepian basis over the rotated domain.
%         - "gl": Use Gauss-Legendre quadrature for the integration (same
%           as in most implementations in the SLEPIAN packages)
%         - "grid": Use a spatial grid to evaluate the Slepian functions
%           and compute the integrals as sums (as described in Bates et al.
%           2017)
%       The Gauss-Legendre method should consume significantly less memory.
%       The default method is "gl".
%   GlNodes (name-value) - Number of Gauss-Legendre quadrature nodes for
%       latitude integration, if the "gl" method is chosen
%       The default number of nodes is 101.
%   GridResFactor (name-value) - Resolution factor for the polar grid used
%       to compute the second localisation matrix (of the polar cap Slepian
%       basis), if the "grid" method is chosen
%       The integral for the localisation matrix will be evaluated over a
%       polar grid with roughly L^2 * GridResFactor points.
%       The default value is 8.
%   ForceNew (name-value) - Force recomputation and overwrite of save files
%       The default value is FALSE.
%   SaveData (name-value) - Save the resulting basis to IFILES path
%       The default value is TRUE.
%   BeQuiet (name-value) - Suppress console output messages
%       The default value is FALSE.
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
%   K - Localisation matrix for the polar cap Slepian basis over the
%       rotated domain
%       This is always computed for the full polar-cap Slepian basis,
%       independent of any truncation applied to G or V.
%       Size: [numPCFuns x numPCFuns], where numPCFuns is the number of
%       polar-cap Slepian functions before truncation.
%
% See also
%   GLMALPHA, GRUNBAUM, KERNELCP
%
% Author
%	2026/03/05, En-Chi Lee (williameclee@arizona.edu)
%
% Last modified
%	2026/03/20, En-Chi Lee (williameclee@arizona.edu)
%     - Saved results to disc for reuse
%     - Made the localisation matrix an output argument
%     - Added the SLEPIAN_ALPHA way of integrating the localisation matrix
%       and made it the default method
%     - Modularised the localisation matrix computation
%	2026/03/06, En-Chi Lee (williameclee@arizona.edu)
%     - Added better guards and log messages for truncation and domain
%       containment
%     - Changed eigenvalues (V) output format
%     - Added truncation and rotb arguments

function [G, V, N, K] = glmalpha_eff(domain, L, truncation, rotb, options)

    arguments (Input)
        domain
        L (1, 1) {mustBeInteger, mustBePositive} = 18
        truncation (1, 1) = NaN
        rotb (1, 1) {mustBeNumericOrLogical} = true
        options.pcapConcThreshold (1, 1) ...
            {mustBeInRange(options.pcapConcThreshold, 0, 1, "exclude-lower")} = 0.3
        options.IntegrationMethod ...
            {mustBeTextScalar, mustBeMember(options.IntegrationMethod, ["grid", "gl"])} = "gl"
        options.GlNodes (1, 1) {mustBePositive, mustBeInteger} = 101
        options.GridResFactor (1, 1) {mustBePositive} = 8
        options.ForceNew (1, 1) logical = false
        options.SaveData (1, 1) logical = true
        options.BeQuiet (1, 1) logical = false
    end

    arguments (Output)
        G (:, :) {mustBeReal, mustBeFinite}
        V (2, :) {mustBeNonnegative}
        N (1, 1) {mustBePositive}
        K (:, :) {mustBeReal, mustBeFinite}
    end

    if isnumeric(truncation)

        if ~isnan(truncation) && (truncation <= 0)
            error("Truncation must be positive, but got %s.", num2str(truncation));
        end

    elseif (isstring(truncation) || ischar(truncation))

        if ~strcmpi(truncation, "N")
            error("Truncation must be either a positive value or the string 'N', but got '%s'.", ...
                truncation);
        end

    else
        error("Truncation must be either a positive value or the string 'N', but got class %s.", ...
            upper(class(truncation)));
    end

    pcapConcThreshold = options.pcapConcThreshold;

    %% Prcoputation check
    dataPath = getoutputfile(domain, L, pcapConcThreshold, rotb, options);

    vars = {'G', 'V', 'N', 'K'};

    if ~options.ForceNew && exist(dataPath, 'file') && ...
            all(ismember(vars, who('-file', dataPath)))
        data = load(dataPath, vars{:});
        G = data.G;
        V = data.V;
        N = data.N;
        K = data.K;

        if ~options.BeQuiet
            fprintf('[SLEPIAN>%s] Loaded <a href="matlab: fprintf(''%s\\n'');open(''%s'')">efficient Slepian projection matrix</a>.\n', ...
                mfilename, dataPath, dataPath);
        end

        % Truncate the basis if asked
        if (isstring(truncation) || ischar(truncation)) && strcmpi(truncation, "N")
            truncation = round(N);
        end

        if ~isnan(truncation)
            truncation = round(truncation); % In case it's a non-integer numeric value

            if truncation <= size(G, 2)
                G = G(:, 1:truncation);
                V = V(:, 1:truncation);
            else
                warning( ...
                    ['Truncation level (%d) is larger than the number of loaded functions (%d). ', ...
                 'No truncation applied.'], ...
                    truncation, size(G, 2));
            end

        end

        return
    end

    %% Main computation
    tic;

    if ~options.BeQuiet
        fprintf( ...
            ['[SLEPIAN>%s] Computing Slepian basis with the efficient formulation, ', ...
         'this may take a while...\n'], ...
            mfilename);
    end

    % Step 1: Find the enclosing polar cap
    [pcapLonlatd, ~] = enclosingCap(domain, "OutputUnit", "degrees");

    % Step 2: Rotate the domain to the North Pole
    pLonlatd = rotateToNPole(domain, pcapLonlatd, "OutputUnit", "degrees");
    % Make sure the domain is actually enclosed
    radiusd = 90 - min(pLonlatd(:, 2));

    % Step 3: compute the Slepian functions for the polar cap
    pcapConcs = []; % The eigenvalues
    pcapGs = {}; % The Slepian coefficients
    pcapMs = []; % The order each function corresponds to

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

        % Step 3c: Evaluate the Slepian functions spatially is absorbed into step 4
    end

    % Sort the Slepian functions by concentration
    [pcapConcs, pcapConcSortId] = sort(pcapConcs, "descend");
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
    switch options.IntegrationMethod
        case "gl"
            K = localisationMatrix(...
                L, pcapGs, pcapMs, radiusd, pLonlatd, options.GlNodes);
        case "grid"
            K = localisationMatrix_grid(...
                L, pcapGs, pcapMs, radiusd, pLonlatd, options.GridResFactor);
        otherwise
            error("slepian:efficientSlepian:invalidIntegrationMethod", ...
                ['Integration method must be either "gl" or "grid", ', ...
                'but got invalid option "%s".'], ...
                options.IntegrationMethod);
    end

    % Step 5: Eigen-decomposition of localisation matrix
    % Get the Slepian functions for the polar cap Slepian functions
    [pSlepG, pSlepConcs] = eig(K);
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

    if ~options.BeQuiet
        t = toc;
        fprintf('[SLEPIAN>%s] Finished computing Slepian basis in %.2f seconds.\n', ...
            mfilename, t);
    end

    N = (L + 1) ^ 2 * spharea(pLonlatd);
    V = [pcapConcs(:), pSlepConcs(:)].'; % Concentrations (eigenvalues)

    if any(V > 1, "all")
        warning("slepian:efficientSlepian:invalidEigenvalues", ...
            '%d eigenvalues are greater than 1 (max: %.3f), which should not happen.', ...
            sum(V > 1, "all"), max(V, [], "all"));
    end

    if options.SaveData
        save(dataPath, '-v7.3', 'G', 'V', 'N', 'K');

        if ~options.BeQuiet
            fprintf('[SLEPIAN>%s] Saved <a href="matlab: fprintf(''%s\\n'');open(''%s'')">efficient Slepian projection matrix</a>.\n', ...
                mfilename, dataPath, dataPath);
        end

    end

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
                ['Truncation level (%d) is larger than the number of computed functions (%d). ', ...
             'No truncation applied.'], ...
                truncation, numFuns);
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

function locMat = ...
        localisationMatrix(L, pcapGs, pcapMs, radiusd, pLonlatd, nGL)
    % Computes the localisation matrix for the polar cap Slepian basis over
    % the rotated domain, using the same method as in KERNELCP

    arguments (Output)
        locMat (:, :) {mustBeReal, mustBeFinite}
    end

    % Gauss-Legendre quadrature interval and nodes over the colatitude range of the polar cap
    [glWeights, glNodes, ~] = gausslegendrecof(nGL, [], ...
        [cosd(radiusd), cosd(0)]);

    % Evaluate the colatitude profiles
    uniqueMs = unique(pcapMs);
    numFuns = length(pcapMs);
    % The Slepian colatitude evaluations at GL nodes [nodes x funs]
    pcapSlepColats = nan(length(glNodes), numFuns);

    for im = 1:length(uniqueMs)
        m = uniqueMs(im);
        idx = (pcapMs == m);
        pcapG_m = pcapGs(idx);

        [Xlm, ~, ~] = xlm(abs(m):L, abs(m), acos(glNodes), 0);

        % Format into [ells x nodes]
        Xlm = reshape(Xlm, [], length(glNodes));

        % Colatitude profile for each retained function, [nodes x funs]
        % Since the zonal component of SHs making up a Slepian function is
        % just the same sine or cosine, we can sum the longitudinal part of
        % the Ylm (i.e. Xlm) and then later multiply by the zonal part
        pcapSlepColat_m = Xlm' * cell2mat(pcapG_m(:)');

        % Append to our collection
        pcapSlepColats(:, idx) = pcapSlepColat_m;
    end

    % Longitudinal integration intervals for the domain at the GL nodes
    zonalIntervals = deg2rad(dphregion(acosd(glNodes), [], pLonlatd));

    locMat = nan(numFuns, numFuns);

    for i = 1:numFuns

        for j = i:numFuns
            m1 = pcapMs(i);
            m2 = pcapMs(j);

            if m1 > 0 && m2 > 0
                zonalIntg = sinsin(acos(glNodes), m1, m2, zonalIntervals);
            elseif m1 <= 0 && m2 <= 0
                zonalIntg = coscos(acos(glNodes), m1, m2, zonalIntervals);
            elseif m1 > 0 && m2 <= 0
                zonalIntg = sincos(acos(glNodes), m1, m2, zonalIntervals);
            else
                zonalIntg = sincos(acos(glNodes), m2, m1, zonalIntervals);
            end

            % Apply normalisation factor for non-zonal orders
            % Ylm = Xlm * sqrt(2 - (m==0)) * trig(m * phi)
            normFactor = sqrt(2 - (m1 == 0)) * sqrt(2 - (m2 == 0));
            zonalIntg = zonalIntg * normFactor;

            locMat(i, j) = sum(glWeights(:)' .* ...
                pcapSlepColats(:, i)' .* pcapSlepColats(:, j)' .* zonalIntg(:)');
            locMat(j, i) = locMat(i, j); % symmetric
        end

    end

end

function locMat = ...
        localisationMatrix_grid(L, pcapGs, pcapMs, radiusd, pLonlatd, res)
    % Computes the localisation matrix for the polar cap Slepian basis over
    % the rotated domain, using the spatial grid method

    % The spatial grid to evaluate the Slepian functions over
    [pgridLond, pgridLatd, ~, pgridWeight, pgridMask] = ...
        polarGridMask(radiusd, pLonlatd, L, resFactor = res);
    pgridWeight = pgridWeight .* pgridMask; % Mask the weights

    % Evaluate the Slepian functions for the polar cap basis over the grid
    uniqueMs = unique(pcapMs);
    numFuns = length(pcapMs);

    pcapSlep = nan(length(pgridLatd), length(pgridLond), numFuns);

    for im = 1:length(uniqueMs)
        m = uniqueMs(im);
        idx = (pcapMs == m);
        pcapG_m = pcapGs(idx);

        Ylm_m = zeros(length(pgridLatd), length(pgridLond), L - abs(m) + 1);

        for l = abs(m):L
            Ylm_m(:, :, l - abs(m) + 1) = ...
                ylm(l, m, deg2rad(90 - pgridLatd), deg2rad(pgridLond));
        end

        id = find(idx);

        for i = 1:length(id)
            f = id(i);
            pcapSlep_i = sum(Ylm_m .* ...
                reshape(pcapG_m{i}, 1, 1, []), 3);
            pcapSlep(:, :, f) = pcapSlep_i;
        end

    end

    % Computes the inner products of the Slepian functions
    locMat = nan(numFuns, numFuns);

    for i = 1:numFuns
        locMat(i, i) = sum(pcapSlep(:, :, i) .^ 2 .* pgridWeight, "all");

        for j = i + 1:numFuns
            locMat(i, j) = sum( ...
                pcapSlep(:, :, i) .* pcapSlep(:, :, j) .* pgridWeight, "all");
            locMat(j, i) = locMat(i, j);
        end

    end

end

function dataPath = getoutputfile(domain, L, pcapConcThreshold, rotb, options)
    % Generates the human-readable string for the saved matrix filename
    if isa(domain, 'char') || isa(domain, 'string')
        domainId = char(domain);
    elseif isa(domain, 'GeoDomain')
        domainId = char(domain.Id);
    elseif iscell(domain)
        domainId = char(domain{1});
    else
        domainId = hash(domain, 'sha1');
    end

    switch options.IntegrationMethod
        case "gl"
            intStr = sprintf('gl_%i', options.GlNodes);
        case "grid"
            intStr = sprintf('grid_%g', options.GridResFactor);
    end

    if ~rotb
        rotbStr = '-norot';
    else
        rotbStr = '';
    end

    outputFile = sprintf('glmalpha_eff-%s-%i-%g-%s%s.mat', ...
        domainId, L, pcapConcThreshold, intStr, rotbStr);

    dataFolder = fullfile(getenv('IFILES'), 'GLMALPHA_EFF');

    if ~exist(dataFolder, 'dir')
        mkdir(dataFolder);
    end

    dataPath = fullfile(dataFolder, outputFile);
end
