function V_max = compute_vmax_prop(D, RPM, J_max, M_tip_max, a)
% COMPUTE_VMAX_PROP  Determine the forward speed limit of a propeller due to
%                    (1) blade‑tip Mach number and (2) maximum advance ratio.
%
%   V_max = compute_vmax_prop(D, RPM, J_max, M_tip_max, a)
%
% INPUTS
%   D         : Propeller diameter [m]
%   RPM       : Rotational speed [rev/min]
%   J_max     : Maximum usable advance ratio (= V / nD) where thrust → 0
%   M_tip_max : Allowable blade‑tip Mach number (e.g. 0.85 for fixed‑pitch)
%   a         : Speed of sound at operating altitude [m/s]
%
% OUTPUT
%   V_max     : The lesser of the Mach‑limited and advance‑ratio‑limited
%               forward speeds [m/s]
%
% METHOD
%   1. Convert RPM to angular velocity (rad/s).
%   2. Compute the rotational component of blade‑tip speed.
%   3. Deduce the translational component permitted before exceeding the
%      Mach limit (Pythagorean addition of velocities).
%   4. Compute the speed corresponding to the maximum advance ratio.
%   5. Return the lower of the two limits as V_max.
% ------------------------------------------------------------------------

%% 1) Convert RPM → angular speed & geometry ------------------------------
n     = RPM / 60;          % Shaft speed [rev/s]
Omega = 2 * pi * n;        % Angular speed [rad/s]
R     = D / 2;             % Prop radius [m]

%% 2) Mach‑number constraint --------------------------------------------
U_rot      = Omega * R;        % Purely rotational tip speed [m/s]
U_tip_max  = M_tip_max * a;    % Absolute tip‑speed limit from Mach cap

% Pythagorean relation: (U_rot)^2 + (V)^2 ≤ (U_tip_max)^2
V_sq = U_tip_max^2 - U_rot^2;  % Remaining margin for axial flight speed
if V_sq < 0
    V_max_Mach = 0;            % Already over‑Mach while static
else
    V_max_Mach = sqrt(V_sq);
end

%% 3) Advance‑ratio (thrust) constraint ----------------------------------
% Advance ratio J = V / (n D).  Beyond J_max, thrust coefficient ≈ 0.
V_max_J = J_max * n * D;

%% 4) Select the governing limit -----------------------------------------
V_max = min(V_max_Mach, V_max_J);

%% 5) User‑friendly printout ---------------------------------------------
fprintf('V_{max} (Mach limit)    = %6.2f m/s  (%4.1f kt, %4.1f mph)\n', ...
        V_max_Mach, V_max_Mach * 1.94384, V_max_Mach * 2.23694);
fprintf('V_{max} (Adv. ratio)    = %6.2f m/s  (%4.1f kt, %4.1f mph)\n', ...
        V_max_J,    V_max_J    * 1.94384, V_max_J    * 2.23694);
fprintf('Final V_{max}           = %6.2f m/s  (%4.1f kt, %4.1f mph)\n', ...
        V_max,      V_max      * 1.94384, V_max      * 2.23694);

end % ========================= END OF FUNCTION ==========================
