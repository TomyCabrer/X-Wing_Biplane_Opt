function [c, w_wing] = sizeWingChord(b, ConstWeight, V_flight, CL_loop, rho, geom)
%SIZEWINGCHORD   Determine wing chord that meets lift & weight targets.
% -------------------------------------------------------------------------
%  [c, w_wing] = sizeWingChord(b, ConstWeight, V_flight, CL_loop, rho, geom)
%
%  Solves for the single wing‑chord length *c* that allows a rectangular
%  (constant‑chord) wing of span *b* to generate enough lift at the chosen
%   flight speed *V_flight* while simultaneously satisfying the
%  current lift‑coefficient schedule *CL_loop*.  The resulting chord is fed
%  into the structural weight model `wing_weight_finder` to obtain the wing
%  weight *w_wing*.
%
%  INPUTS
%    b            – wingspan                                                  [m]
%    ConstWeight  – aircraft weight that must be supported (N, not kg!)       [N]
%    V_flight     – flight speed for lift sizing            [m/s]
%    CL_loop      – lift coefficient required at the sizing point             [-]
%    rho          – air density                                               [kg/m³]
%    geom         – struct of geometric & material parameters used by
%                   wing_weight_finder (see main script for field list)
%
%  OUTPUTS
%    c            – chord length that balances lift & weight                  [m]
%    w_wing       – structural wing weight estimated by wing_weight_finder    [kg]
%
%  METHOD
%    • Assume a rectangular wing (constant chord) as first approximation.
%    • Formulate a residual function whose root corresponds to the chord
%      length that produces the required lift:  L(c) – Weight = 0.
%    • Bracket the root around an initial guess by scaling the chord by an
%      order of magnitude either side; if no sign change is found, the
%      function exits with NaN results.
%    • Solve the scalar nonlinear equation with `fzero`.  A tight tolerance
%      (`TolX = 1e-10`) avoids numerical noise propagating downstream.
%    • Feed the converged chord into the structural model to compute wing
%      mass.
% -------------------------------------------------------------------------
%  NOTE
%    The variable `ok` is retained from the legacy implementation but is
%    *not* returned; it merely serves local error‑handling logic.
% -------------------------------------------------------------------------

    % === 1) Initial guess -------------------------------------------------
    AR0 = 6;                 % Choose a moderate aspect‑ratio starting point
    c0  = b / AR0;           % First‑guess chord (m)

    % === 2) Define residual function -------------------------------------
    % The nested handle calls the helper `chordResidual`, which encapsulates
    % lift balance + any geometry‑dependent effects captured in `geom`.
    chordFun = @(c) chordResidual(c, b, ConstWeight, V_flight, CL_loop, rho, geom);

    % === 3) Bracket a sign change around c0 ------------------------------
    % Generate a lower & upper bound one order of magnitude apart so that
    % `fzero` is guaranteed a valid sign change interval.
    br = [0.2*c0, 5*c0];     % [m] bracket vector
    F1 = chordFun(br(1));
    F2 = chordFun(br(2));

    % Abort early if lift residual does not change sign or returns NaN/Inf.
    if any(~isfinite([F1, F2])) || sign(F1) * sign(F2) >= 0
        c = NaN;  w_wing = NaN;  ok = 0;  return
    end

    % === 4) Root‑find inside the bracket ---------------------------------
    try
        opts        = optimset('TolX', 1e-10, 'Display', 'off');
        [c, ~, eflag] = fzero(chordFun, br, opts);
    catch           % Defensive on unexpected optimiser errors
        eflag = -1;
    end

    % Validate solver outcome before proceeding
    if eflag <= 0 || ~isfinite(c) || c <= 0
        c = NaN;  w_wing = NaN;  ok = 0;  return
    end

    % === 5) Structural weight estimate -----------------------------------
    w_wing = wing_weight_finder(geom.top, geom.bottom, b, c, geom.length_x, ...
        geom.xq, geom.dx, geom.Cross_section_N, geom.Spar_length, ...
        geom.aero_thickness, geom.skin_thickness, ...
        geom.rho_spar, geom.rho_skin, geom.rho_rib, geom.thetaDeg);

    ok = 1;   % (not returned)
end
