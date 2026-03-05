%% GRUNBAUM - Calculates the matrix the way Grunbaum et al. (1982) propose
% Eigenfunctions of the fixed-order SINGLE POLAR CAP concentration problem.
% Orders the eigenfunctions in decreasing order.
%
% Syntax
%   [E, Vg, th, C, T, V] = GRUNBAUM(TH, L, m, nth, grd)
%
% INPUT:
% TH          Angular extent of the spherical cap, in degrees
% L           Bandwidth (maximum angular degree)
% m           Angular order of the required data window, -l<m<l
% nth         Number of points sampling between 0 and pi [default: 720]
%             If nth=0, E will be zero and no sign-correction performed
% grd         1 Colatitudes only; returns matrix E [default]
%             2 Colatitude/Longitude; returns cell E
% xver        0 Do nothing else [default: 0]
%             1 Excessive verification
%
% OUTPUT:
%
% E           The spatial eigenfunctions, not multiplied by sin or cos(mphi)
% Vg          The Grunbaum eigenvalues
% th          The colatitudes at which the functions are evaluated
% C           The real spherical harmonic expansion coefficients
% T           The tridiagonal matrix
% V           The eigenvalues given by integration over the patch
%
% EXAMPLE
%   Compare Grunbaum's faster approach with our numerical slower one
%   [E,Vg,th,C,T,V]=grunbaum(20,18,0,180);
%   [E2,V2,N,th,C2,~,~,~,~,~,K]=sdwcap(20,18,0,180);
%   difer(abs(T*K-K*T))
%
% SEE ALSO
%   SDWCAP, BOXCAP
%
% Last modified
%	2026/03/04, En-Chi Lee (williameclee@arizona.edu)
%     - Small refactors, should not change any behaviour
%   2012/03/21, FJ Simons (fjsimons@alum.mit.edu)

function [E, Vg, th, C, T, V] = grunbaum_new(radius, L, m, nth, grd, xver)

    arguments (Input)
        radius (1, 1) {mustBeNumeric, mustBePositive} = 40
        L (1, 1) {mustBeNumeric, mustBeInteger} = 18
        m (1, 1) {mustBeNumeric, mustBeInteger} = 0
        nth (1, 1) {mustBeNumeric} = 720
        grd (1, 1) {mustBeNumeric} = 0
        xver (1, 1) {mustBeNumeric} = 0
    end

    if m > L || m < -L
        error('Order cannot exceed degree (i.e. must have %d <= m <= %d), but got m = %d', -L, L, m)
    end

    mor = m;
    m = abs(m);

    n = L - m + 1;
    b = cosd(radius);

    % This is Grunbaum's original result, only for m>=0
    e = 1:n;
    alpha = -(e + m) .* (e + m - 1);
    gamma = sqrt((e + m) .^ 2 - m ^ 2) .* (((e + m) .^ 2 - (L + 1) ^ 2) ./ sqrt(4 * (e + m) .^ 2 - 1));
    T = tridiag(gamma(1:end - 1), alpha * b, gamma(1:end - 1));

    if xver == 1
        % This result should be completely equivalent
        e = m:L;
        Gll = -e .* (e + 1);
        Gl1 = (e .* (e + 2) - L * (L + 2)) .* sqrt(((e + 1) .^ 2 - m ^ 2) ./ (2 * e + 1) ./ (2 * e + 3));
        TT = tridiag(Gl1(1:end - 1), Gll * b, Gl1(1:end - 1), n);
        difer(sum(T(:) - TT(:)))
    end

    % Diagonalise the matrix
    [C, Vg] = eig(T);

    % Turn this into a vector so orthocheck knows the concentration factors
    % still need to be calculated.
    Vg = diag(Vg);

    % Check normalization and calculate the eigenvalues
    % from a straightforward GL integration
    [~, ~, ~, V] = orthocheck(C, Vg, deg2rad(radius), m);

    if nth == 0
        E = 0;
        th = 0;
        return
    end

    % Compute spatial functions, colatitudinal part only
    % Zonal functions only
    if m == 0
        % Make spatial functions
        % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
        % (5.12) already included!
        [E, th] = pl2th(C, nth, 1);
        th = th * 180 / pi;
        nlon = 2 * nth - 1;
    else
        % This is SDW (2005) equation (5.10) combined with the sqrt(2-dom) of
        % (5.12) already included!
        [E, nlon, ~] = plm2th(C, nth, m, 1);
        th = linspace(0, 180, size(E, 1));
    end

    % Make E start with a positive lobe and ajust C too
    % Take the first NONZERO sample! Not just a numbered sample!
    % This was the source of a very nasty bug
    for index = 1:size(E, 2)
        C(:, index) = C(:, index) * sign(indeks(E(~ ~E(:, index), index), 1));
        E(:, index) = E(:, index) * sign(indeks(E(~ ~E(:, index), index), 1));
    end

    if grd == 2
        % Output on full grid; watch the sign of m
        % Negative order is cosine taper
        if mor <= 0
            EE = E; clear E

            for index = 1:size(EE, 2)
                E{index} = EE(:, index) * cos(m * linspace(0, 2 * pi, nlon));
            end

        end

        % Positive order is sine taper
        if mor > 0
            EE = E; clear E

            for index = 1:size(EE, 2)
                E{index} = EE(:, index) * sin(m * linspace(0, 2 * pi, nlon));
            end

        end

    end

end
