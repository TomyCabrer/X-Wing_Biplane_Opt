function F = chordResidual(c, b, ConstWeight, V_flight, CL_loop, rho, geom)
%CHORDRESIDUAL  Residual function for wing-chord root‑finding.
% -------------------------------------------------------------------------
%   F = chordResidual(c, b, ConstWeight, V_flight, CL_loop, rho, geom)
%
%   Evaluates the scalar residual used by `sizeWingChord` when solving for
%   the constant wing‑chord *c* that balances lift with aircraft weight.
%   The residual is simply the difference between the chord required by the
%   lift equation (*c_req*) and the current trial value (*c*).  A root of
%   this function therefore corresponds to a chord that yields exactly the
%   lift needed to support the total aircraft weight at the minimum
%   sizing speed.
%
%   INPUTS
%     c            – trial chord length being tested                          [m]
%     b            – fixed wingspan                                           [m]
%     ConstWeight  – constant aircraft weight excluding the wing              [N]
%     V_flight     – sizing flight speed                                      [m/s]
%     CL_loop      – lift coefficient required at V_flight                    [-]
%     rho          – air density                                              [kg/m³]
%     geom         – geometry & material struct forwarded to
%                    `wing_weight_finder` (see caller for fields)
%
%   OUTPUT
%     F            – scalar residual  F(c) = c_req(c) – c                     [m]
%                    Root (F = 0) implies lift balance.
%
%   METHOD
%     1) Compute the wing structural weight for the current chord via
%        `wing_weight_finder`.
%     2) Add that to the fixed aircraft weight to obtain the total weight
%        that must be supported.
%     3) Use the steady‑level flight lift equation to back‑solve the chord
%        required to generate that lift at V_flight.
%     4) Return the difference between the required chord and the guess.
% -------------------------------------------------------------------------
%   NOTE
%     If this routine is called repeatedly inside a solver (e.g. `fzero`),
%     runtime can be dominated by the structural weight call.  Once the
%     solver converges successfully, caching the final `w_wing` and using
%     it as the initial guess in subsequent design loops may accelerate
%     higher‑level optimisation passes.
% -------------------------------------------------------------------------

    % === 1) Structural weight for current chord --------------------------
    w_wing = wing_weight_finder(geom.top, geom.bottom, b, c, geom.length_x, ...
        geom.xq, geom.dx, geom.Cross_section_N, geom.Spar_length, ...
        geom.aero_thickness, geom.skin_thickness, ...
        geom.rho_spar, geom.rho_skin, geom.rho_rib, geom.thetaDeg);

    % === 2) Total aircraft weight ----------------------------------------
    Wtot = ConstWeight + w_wing;     % [N] includes wing weight just computed

    % === 3) Chord required by lift equation ------------------------------
    %   Lift equation rearranged for chord (assuming rectangular planform):
    %       L = 0.5 * rho * V^2 * CL * S ,  with  S = b * c
    %       ⇒ c_req = Wtot / (rho * V^2 * CL * b)
    c_req = Wtot ./ (rho .* V_flight.^2 .* CL_loop .* b);

    % === 4) Residual ------------------------------------------------------
    %   Positive F means trial chord is too small (requires larger chord).
    %   Negative F means trial chord is too large.
    F = c_req - c;
end
