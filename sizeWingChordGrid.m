function [C, M] = sizeWingChordGrid(b, ConstWeight, Vmin_lift, CL_loop, rho, geom)
% SIZEWINGCHORDGRID  Vectorised wrapper around sizeWingChord() that computes
%                    the mean aerodynamic chord (C) and wing mass (M) for a
%                    lattice of span \times weight combinations.
%
%   [C, M] = sizeWingChordGrid(b_vec, ConstWeight_vec, Vmin_lift, ...
%                              CL_loop, rho, geom)
%
% INPUTS
%   b_vec        : Column vector (c×1) of candidate half-spans [m]. Each
%                  element represents a distinct wing planform width at the
%                  current design point.
%   ConstWeight  : Row vector (1×z) of aircraft baseline weights [N] that
%                  exclude wing structure mass (thus, a sweep over payload
%                  or battery configurations).
%   Vmin_lift    : Minimum flight speed providing sufficient lift [m/s] for
%                  the current CL_loop value. Passed directly to the core
%                  sizing routine.
%   CL_loop      : Aerodynamic lift coefficient under evaluation (scalar).
%   rho          : Ambient air density [kg/m³].
%   geom         : Auxiliary geometry struct passed through to
%                  sizeWingChord() (e.g., airfoil coordinates, spar layout).
%
% OUTPUTS
%   C            : (c×z) matrix of mean chord lengths [m] – each row
%                  corresponds to a span entry in b_vec, each column to a
%                  ConstWeight entry.
%   M            : (c×z) matrix of wing structural mass [kg] matching C.
%
% METHOD
%   The function simply iterates over every (span, weight) pair and calls
%   the more detailed sizeWingChord() routine that solves the local
%   structural/aerodynamic sizing problem. Although nested loops are used,
%   the dimensionality is modest (O(10²)), so explicit looping keeps the
%   code readable without incurring noticeable runtime penalties.
% -------------------------------------------------------------------------

% Ensure consistent orientation of input vectors -------------------------
b            = b(:);        % Force column vector for spans (c×1)
ConstWeight  = ConstWeight(:).';   % Force row vector for weights (1×z)

c = numel(b);              % Number of span candidates
z = numel(ConstWeight);    % Number of weight candidates

% Pre-allocate result matrices for speed ---------------------------------
C = zeros(c, z);           % Mean chord matrix [m]
M = zeros(c, z);           % Wing mass  matrix [kg]

% Nested sweep over design grid ------------------------------------------
for i = 1:c                % Loop over spans
    for j = 1:z            % Loop over weights
        [C(i,j), M(i,j)] = sizeWingChord( ...
            b(i), ConstWeight(j), Vmin_lift, CL_loop, rho, geom);
    end
end

end % ========================= END OF FUNCTION ==========================
