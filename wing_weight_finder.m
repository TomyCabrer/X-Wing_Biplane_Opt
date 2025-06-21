function wing_weight = wing_weight_finder(top, bottom, b, c, length_x, ...
                                         xq, dx, Cross_section_N, Spar_length, ...
                                         aero_thickness, skin_thickness, ...
                                         rho_spar, rho_skin, rho_rib, thetaDeg)
% WING_WEIGHT_FINDER  Estimate the structural weight (N) for a semi‑monocoque
%                     wing box composed of spars, ribs (cross‑sections) and
%                     skin panels.
%
%   wing_weight = wing_weight_finder(top, bottom, b, c, length_x, xq, dx, ...
%                                    Cross_section_N, Spar_length, ...
%                                    aero_thickness, skin_thickness, ...
%                                    rho_spar, rho_skin, rho_rib, thetaDeg)
%
% INPUTS
%   top, bottom     : Nx2 arrays containing the airfoil upper and lower
%                     surfaces sampled along the unit‑chord (x, y).
%   b               : Wing semi‑span [m].
%   c               : Mean aerodynamic chord [m].
%   length_x        : Span portion inboard of the dihedral break [m].
%   xq              : Vector of x‑coordinates (unit‑chord) used for surface
%                     integration.
%   dx              : Scalar spacing between consecutive xq samples.
%   Cross_section_N : Number of full rib stations per semi‑wing.
%   Spar_length     : Depth of main spar web/boom [m] (assumed square).
%   aero_thickness  : Fractional airfoil thickness (t/c).
%   skin_thickness  : Panel thickness [m].
%   rho_*           : Densities for spar, skin, and rib materials [kg/m³].
%   thetaDeg        : Wing dihedral angle [deg].
%
% OUTPUT
%   wing_weight     : Total structural weight for one semi‑wing [N].
%
% NOTE
%   * The method assumes a simple prismatic spar (square cross‑section),
%     identical ribs at evenly‑spaced stations, and a lofted skin surface.
%   * All three major contributors are multiplied by two to account for
%     symmetric left/right semi‑wings.
% -------------------------------------------------------------------------

%% 1) MAIN SPAR MASS -------------------------------------------------------
% Square spar cross‑sectional area [m²]
sparA = Spar_length^2;

% Spar volume: two spars run along both the inboard and outboard portions
% of the wing (factor 2), scaled for dihedral on the inner panel.
sparV = 2 * sparA .* (length_x ./ cosd(thetaDeg) + (b - length_x));

% Mass [kg]
sparM = sparV .* rho_spar;

%% 2) RIB (AIRFOIL CROSS‑SECTION) MASS ------------------------------------
% Integrate area enclosed between upper and lower surfaces on unit chord
crossA = 0;
for i = 2:length(xq)
    crossA = crossA + dx .* (top(i,2) - bottom(i,2));
end

% Volume of ribs: 4 to cover both left/right wings and front/back halves
crossV = 4 .* Cross_section_N .* crossA .* c.^2 .* aero_thickness;

% Mass [kg]
crossM = crossV .* rho_rib;

%% 3) SKIN MASS -----------------------------------------------------------
% Integrate surface length of both upper and lower skins on unit chord
s = 0;
for i = 2:length(xq)
    s = s + sqrt(dx.^2 + (top(i,2)    - top(i-1,2)).^2);   % Upper surface
    s = s + sqrt(dx.^2 + (bottom(i,2) - bottom(i-1,2)).^2);% Lower surface
end
s = s .* c;  % Dimensionalise by chord length

% Total skin area covering both panels (factor of 2 for symmetric wings)
skinA = 2 .* s .* (length_x ./ cosd(thetaDeg) + (b - length_x));

% Volume & mass of skin
skinV = skinA .* skin_thickness;
skinM = skinV .* rho_skin;

%% 4) TOTAL WING WEIGHT ----------------------------------------------------
wing_weight = (sparM + crossM + skinM) .* 9.81;   % Convert kg → N

end % ========================= END OF FUNCTION ==========================
