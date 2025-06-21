function [CL_data,bestCL] = OptFun(Power_W, Thrust_N, Throttle, n_motors, ...
                                    CL, V_input, ConstWeight, battery_mass, ...
                                    b, battery_capacity, rho, cmin, h, factor, ...
                                    kM, dF, mu, ThrustMT, AoA, thetaDeg, geom, ...
                                    max_speed, length_x, Payload_mass)
% OPTFUN  Evaluates mission performance over a grid of lift coefficients and
%         flight velocities and returns the Pareto-optimal combination.
%
% USAGE:
%   [CL_data, bestCL] = OptFun(Power_W, Thrust_N, Throttle, n_motors, ...
%                              CL, V_input, ConstWeight, battery_mass, ...
%                              b, battery_capacity, rho, cmin, h, factor, ...
%                              kM, dF, mu, ThrustMT, AoA, thetaDeg, geom, ...
%                              max_speed, length_x, Payload_mass)
%
% INPUTS
%   Power_W         : Vector of motor shaft power [W] @ discrete thrust pts
%   Thrust_N        : Corresponding thrust vector per motor [N]
%   Throttle        : ESC throttle map (fraction) for each thrust point
%   n_motors        : Number of propulsion units installed
%   CL              : Grid of candidate lift-coefficients to sweep through
%   V_input         : Grid of candidate flight speeds to sweep [m/s]
%   ConstWeight     : Baseline aircraft weight excl. wings & battery [N]
%   battery_mass    : Vector of discrete battery masses considered [kg]
%   b               : Vector of wingspans matching CL grid [m]
%   battery_capacity: Energy content for each battery mass element [J]
%   rho             : Air-density at design altitude [kg/m³]
%   cmin            : Minimum allowable mean chord [m] (structural limit)
%   h               : Fuselage height for wetted area calc [m]
%   factor          : Oswald-efficiency factor tuning constant (empirical)
%   kM              : Additional span-wise efficiency modifier
%   dF              : Fuselage diameter [m]
%   mu              : Dynamic viscosity of air [Pa·s]
%   ThrustMT        : Maximum thrust available from all motors [N]
%   AoA             : Angle-of-attack grid corresponding to CL [deg]
%   thetaDeg        : Wing dihedral angle [deg]
%   geom            : Additional geometry struct for chord sizing routine
%   max_speed       : Never-exceed air-speed [m/s]
%   length_x        : Wing semi-span over the fuselage [m]
%   Payload_mass    : Fixed payload carried [kg]
%
% OUTPUTS
%   CL_data : Table of non-dominated solutions (rows) with detailed metrics
%   bestCL  : CL value that maximises mission range across full search-space
%
% DESCRIPTION
%   The routine forms a 2-D sweep across lift coefficient (CL) and flight
%   velocity (V_input). For each (CL, V) point it synthesises a matching
%   wing planform, evaluates aerodynamic coefficients (CD0, k, e, CD),
%   calculates thrust & power requirements against the motor map and
%   determines range and endurance from the chosen battery. Infeasible
%   combinations (e.g., exceeding thrust capability, structural limits)
%   are masked with NaNs. Finally, rows that yield the maximum range per
%   CL are extracted and returned as Pareto candidates.
% -----------------------------------------------------------------------

% --------------------------- PRE-ALLOCATE -------------------------------
Range              = zeros(length(CL), length(V_input)); % [m]
time_flight        = Range;    % [min]
Drag               = Range;    % [N]
AspectRatio        = Range;    % [-]
Span               = Range;    % [m]
MassBattery        = Range;    % [kg]
Power_while_flying = Range;    % [W]
CDt                = Range;    % [-]
V_flight           = Range;    % [m/s]
VL                 = Range;    % [m/s] – lift-off speed
TotalMass          = Range;    % [kg]
StoreRange         = 0;        % Keeps record of current best range surface

m     = 1;                               % Row index for CL loop
max_t = max(Thrust_N) * n_motors;        % Absolute thrust ceiling [N]

% ========================= MAIN DESIGN LOOP ============================
for CL_loop = CL'                   % Sweep through candidate lift coeffs
    n = 1;                          % Column index for velocity loop

    for Vmin_lift = V_input         % Sweep through velocities

        % ----------------- Wing sizing for current (CL,V) --------------
        [c, w_wing] = sizeWingChordGrid( ... % Mean chord & wing-mass
            b, ConstWeight, Vmin_lift, CL_loop, rho, geom);

        FinalWeight = ConstWeight + w_wing;   % Updated take-off weight [N]
        c(c < cmin) = NaN;                    % Enforce min chord limit

        AR = b ./ c;                          % Aspect ratio
        V  = Vmin_lift;                       % True flight speed [m/s]

        % ------ Planform area incl. dihedral correction ----------------
        Sref = 2 .* c .* length_x ./ cosd(thetaDeg) ...
             + 2 .* (b - length_x) .* c;

        % -------- Wetted area (fuselage + wings) -----------------------
        Swet = (1.977 + 0.52 * 0.14) .* Sref ...   % Flat plate coeff
             + pi * dF * h ...                    % Fuselage side
             + 0.25 * pi * dF^2 * 2;              % Nose + tail caps

        % --------------------- Parasitic drag --------------------------
        eTheo = (1 + factor * AR) .^ (-1);         % Theoretical Oswald-e
        kF    = 1 - 2 * (dF ./ b).^2;              % Fuselage-span eff
        Re    = (Swet ./ b) .* V .* rho / mu;      % Reynolds number
        Cf    = 0.00258 ...                        % Spalding-Chi drag eqn
              + 0.00102 .* exp(-6.28e-9  .* Re) ...
              + 0.00295 .* exp(-2.01e-8 .* Re);
        CD0   = Cf .* (Swet ./ (2 .* b .* c));     % Zero-lift drag coeff

        % -------- Induced drag & total drag coeff ----------------------
        e  = kM ./ (1 ./ (eTheo .* kF) + 0.38 .* CD0 .* pi .* AR);
        k  = 1 ./ e;                               % Induced drag factor
        CD = CD0 + (CL_loop.^2 .* k) ./ (AR .* pi);

        % ------------------- Thrust & Power req. -----------------------
        D = 0.5 .* rho .* V.^2 .* CD .* Sref;      % Drag force [N]
        D(D > ThrustMT) = NaN;                     % Cap at motor limit

        T = D ./ (1 - V ./ max_speed);             % Thrust inc. inflow
        T(T <= 0 | T > max_t) = NaN;               % Eliminate invalid pts

        % Interpolate electrical power @ required thrust
        power_flight = interp1(Thrust_N .* n_motors, ...
                               Power_W   .* n_motors, ...
                               T, 'pchip');

        % ------------------- Endurance & Range ------------------------
        tf = battery_capacity ./ power_flight;   % [s] endurance
        R  = tf .* V;                            % [m] range

        % Identify best-range operating point within current matrix cell
        [~, idx] = max(R(:));                    % Linear index
        [row, col] = ind2sub(size(R), idx);

        if all(isnan(R), 'all')                  % Infeasible design
            % Store NaNs to preserve matrix shape
            V_flight(m,n)           = NaN;
            time_flight(m,n)        = NaN;
            Range(m,n)              = NaN;
            Drag(m,n)               = NaN;
            AspectRatio(m,n)        = NaN;
            MassBattery(m,n)        = NaN;
            Power_while_flying(m,n) = NaN;
            CDt(m,n)                = NaN;
            VL(m,n)                 = NaN;
            TotalMass(m,n)          = NaN;
            Span(m,n)               = b(row);          % Store span for ref
        else                                    % Feasible point found
            % Extract scalar metrics at arg-max position
            V_flight(m,n)           = V(1);
            time_flight(m,n)        = tf(idx) / 60;      % [min]
            Range(m,n)              = R(idx);
            Drag(m,n)               = D(idx);
            AspectRatio(m,n)        = AR(idx);
            Span(m,n)               = b(row);
            MassBattery(m,n)        = battery_mass(col);
            Power_while_flying(m,n) = power_flight(idx);
            CDt(m,n)                = CD(idx);
            VL(m,n)                 = Vmin_lift;
            TotalMass(m,n)          = FinalWeight(idx);

            % Update global best range map for later contour plotting
            if R(idx) > max(StoreRange(:))
                StoreRange = R;
            end
        end
        n = n + 1;      % Advance velocity column index
    end
    m = m + 1;          % Advance CL row index
end

% =================== POST-PROCESS PARETO CURVE ==========================
R_flat = Range;          % Copy for convenience
R_flat(isnan(R_flat)) = -inf;                 % Replace NaNs → -inf
max_values = max(R_flat, [], 2);              % Best range per CL row
max_values(~isfinite(max_values)) = [];       % Remove rows with all NaNs

% Locate indices of best-range points
pos = arrayfun(@(mv) find(Range == mv, 1, 'first'), max_values); % lin idx
[row_indices, ~] = ind2sub(size(Range), pos);

CL_mat=repmat((CL),1,size(Range,2));
Range_mat=Range;
V_flight_mat=V_flight;

% Extract data corresponding to Pareto front
V_flight = V_flight(pos);
Range    = Range(pos);
time_flight = time_flight(pos);
Drag     = Drag(pos);
AspectRatio = AspectRatio(pos);
Span     = Span(pos);
MassBattery = MassBattery(pos);
Power_while_flying = Power_while_flying(pos);
CDt      = CDt(pos);
TotalMass = TotalMass(pos) / 9.81;           % Convert N → kg

bestCL = CL(Range == max(Range));            % Absolute best CL

% ---------------------- Build output table -----------------------------
Max_Payload = (n_motors * max(Thrust_N)) / 9.81 - TotalMass;  % [kg]

CL         = CL(row_indices);                % Pareto CL list
AoA        = AoA(row_indices);               % Matching AoA list
Chord      = Span ./ AspectRatio;            % Mean aerodynamic chord
Throttle_flight = interp1(Thrust_N .* n_motors, Throttle, Drag, ...
                          'pchip', 'extrap');

CL_data = table(CL, AoA, V_flight, Range, time_flight, Drag, ...
                AspectRatio, Span, Chord, MassBattery, TotalMass, ...
                Max_Payload, Power_while_flying, Throttle_flight, CDt);

% ----------------------- Summary banner -------------------------------
MR = find(Range == max(Range));
fprintf(['\nBest range: %.2f km at %.2f m/s during a %.2f-min flight\n' ...
         'CL = %.4f (AoA %.2f°), take-off mass = %.2f kg, payload margin = %.2f kg\n' ...
         'Battery = %.2f kg, span = %.2f m, chord = %.2f m\n'], ...
        Range(MR)/1000, V_flight(MR), time_flight(MR), ...
        CL(MR), AoA(MR), TotalMass(MR) - Payload_mass, Max_Payload(MR), ...
        MassBattery(MR), Span(MR), Chord(MR));

% -------------------- Visualisation (optional) -------------------------
figure; hold on; grid on;
contourf(V_flight_mat,CL_mat,Range_mat/1000)          % Filled contour of range
cb = colorbar;
cb.Label.String = 'Range (km)';
plot(V_flight, CL, 'r', 'LineWidth', 3);    % Pareto front overlay

xlabel('Flight Velocity (m/s)');
ylabel('Lift Coefficient, C_L');
title('Mission Range across C_L – V Envelope');

set(gca, 'LineWidth', 1.5, 'FontSize', 18);

end % =========================== END OF OPTFUN ==========================
