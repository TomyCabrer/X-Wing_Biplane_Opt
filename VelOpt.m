function [velocity_data,diff_velocity_data]=VelOpt(Power_W,Thrust_N,Throttle,n_motors,...
    bestCL,V_input,ConstWeight,battery_mass,b,battery_capacity,rho,cmin,...
    h,factor,kM,dF,mu,ThrustMT,RPM,thetaDeg,geom,max_speed,length_x)
%VELOPT  Size–Performance sweep over minimum flight velocity.
% =========================================================================
%   [VEL, DIFF] = VELOPT(P, T, Thrtl, Nmot, bestCL, Vmin, W, Mbatt, span, Cbatt, ...)
%
%   **Purpose**
%     Mimics the original `battery.m` but sweeps a *vector* of
%     flight velocities, returning the optimum range / endurance and the
%     associated geometry & power settings for each point.  The routine
%     therefore helps choose the best design point (wing span, battery mass
%     etc.) for the Hermes e‑VTOL concept.
%
%   **Inputs (excerpt)**
%     Power_W          [Nx1]  – Motor electrical power map               (W)
%     Thrust_N         [Nx1]  – Corresponding static thrust map          (N)
%     Throttle         [Nx1]  – Normalised throttle (0–1) for the map    (-)
%     n_motors         scalar – Number of propulsive units               (-)
%     bestCL           scalar – Target lift coefficient in cruise        (-)
%     V_input          [1xM]  – Vector of . flight speeds to sweep    (m/s)
%     ConstWeight      scalar – Aircraft weight excl. wing               (N)
%     battery_mass     [1xK]  – Battery mass design vector               (kg)
%     b                scalar – Initial wing span guess (will adapt)    (m)
%     battery_capacity [1xK]  – Battery energy for each mass            (J)
%     rho, cmin, h, ...         Physical & aero constants               (var)
%
%   **Outputs**
%     velocity_data        – Table of optimum metrics versus V_flight
%     diff_velocity_data   – Table for differential‑thrust post‑analysis
%
%   **Method Overview**
%     1) Loop over each candidate velocity `V`.
%     2) Solve for wing chord `c` such that lift = weight using
%        `sizeWingChordGrid`, respecting a minimum chord constraint.
%     3) Compute wetted area and drag breakdown; constrain to available
%        thrust and propulsive speed limit.
%     4) Interpolate motor map to find power required and thus endurance
%        & range for each battery option; store the *best* range point.
%     5) After the V‑sweep, build summary tables and generate diagnostic
%        plots (range/time, AR vs drag, battery‑mass contours, etc.).
% =========================================================================

% -------------------------------------------------------------------------
% 0) Preparations & pre‑allocation
% -------------------------------------------------------------------------
% ---------------------------- 0. House‑keeping ---------------------------
CL = bestCL;                                 % Target lift coefficient

% Pre‑allocate vectors (same length as Vmin_input)
Range               = zeros(size(V_input))';
time_flight         = Range;  % Endurance (min)
Drag                = Range;  % Drag (N)
AspectRatio         = Range;  % Wing aspect ratio
Span                = Range;  % Wing span (m)
MassBattery         = Range;  % Selected battery mass (kg)
Power_while_flying  = Range;  % Electrical power in cruise (W)
CDt                 = Range;  % Total drag coefficient
V_flight            = Range;  % Cruise speed (m/s)
CD_zero             = Range;  % Zero‑lift CD
K                   = Range;  % Induced‑drag factor 1/eπAR

n = 1;                       % Loop counter for velocity sweep
StoreRange = 0;              % Matrices retained for contour plot
StoreAR    = 0;

max_t = max(Thrust_N)*n_motors;   % Static thrust ceiling

% ----------------------- 1. Loop over candidate Vmin ---------------------
for i=V_input
    
    V = i;                          % Current minimum‑lift velocity (m/s)

    % 1.1  Wing chord sizing so L = W at this Vmin
    [c, w_wing] = sizeWingChordGrid(b, ConstWeight, V, CL, rho, geom);
    FinalWeight = ConstWeight + w_wing; %#ok<NASGU>
    c(c < cmin) = NaN;               % Enforce structural chord limit
    AR = b ./ c;                     % Aspect ratio

    % 1.2  Wetted area (Raymer wing + simple cylindrical fuselage)
    Sref = 2.*c.*length_x./cos(thetaDeg*pi/180) + 2*(b-length_x).*c;
    Swet = (1.977 + 0.52*0.14).*Sref + pi*dF*h + 0.5*pi*dF.^2;

    % 1.3  Zero‑lift drag CD0 via flat‑plate skin‑friction correlation
    eTheo = (1 + factor*AR).^(-1);          % Theoretical Oswald efficiency
    kF    = 1 - 2*(dF./b).^2;               % Fuselage interference factor
    Re    = ((Swet./b).*V)*rho/mu;          % Reynolds number
    Cf    = 0.00258 + 0.00102*exp(-6.28e-9*Re) + 0.00295*exp(-2.01e-8*Re);
    CD0   = Cf.*(Swet./(2.*b.*c));          % Zero‑lift drag coefficient

    % 1.4  Induced drag and total CD
    e  = kM ./ (1./(eTheo.*kF) + 0.38.*CD0.*pi.*AR);  % Oswald incl. CD0 term
    k  = 1./e;                                        % Induced factor
    CD = CD0 + (CL^2 .* k)./(AR.*pi);                 % Total drag coefficient

    % 1.5  Drag force and thrust feasibility (with tip‑speed correction)
    D = 0.5*rho.*V.^2.*CD.*Sref;      % Aerodynamic drag (N)
    D(D > ThrustMT) = NaN;            % Above throttle‑limited thrust?
    T = D ./ (1 - V./max_speed);      % Disk‑loading / tip‑Mach correction
    T(T <= 0 | T > max_t) = NaN;      % Outside physical thrust envelope

    % 1.6  Power & mission (sweep battery masses)
    power_flight = interp1(Thrust_N*n_motors, Power_W*n_motors, T, "pchip");
    tf = battery_capacity ./ power_flight;  % Endurance (s)
    R  = tf .* V;                           % Range (m)

    % 1.7  Pick battery giving max range at this Vmin
    [~, idx] = max(R(:));
    [row, col] = ind2sub(size(R), idx);

    if all(isnan(R), 'all')
        % ── No feasible design point: fill NaNs
        V_flight(n)           = NaN;
        time_flight(n)        = NaN;
        Range(n)              = NaN;
        Drag(n)               = NaN;
        AspectRatio(n)        = NaN;
        Span(n)               = NaN;
        MassBattery(n)        = NaN;
        Power_while_flying(n) = NaN;
        CDt(n)                = NaN;
        CD_zero(n)            = NaN;
        K(n)                  = NaN;
    else
        % ── Store best‑range metrics for this Vmin
        V_flight(n)           = V;
        time_flight(n)        = tf(idx)/60;             % min
        Range(n)              = R(idx);
        Drag(n)               = D(idx);
        AspectRatio(n)        = AR(idx);
        MassBattery(n)        = battery_mass(col);
        Power_while_flying(n) = power_flight(idx);
        CDt(n)                = CD(idx);
        CD_zero(n)            = CD0(idx);
        K(n)                  = k(idx);
        Span(n, 1)            = b(row);

                % Save full matrices if this is the global best range so far
        if R(idx) > max(StoreRange(:))
            StoreRange = R;  StoreAR = AR;
        end
    end

    n = n + 1;   % Next Vmin
end

% ------------------------ 2. Assemble results table ----------------------
Chord = Span ./ AspectRatio;  % Mean aerodynamic chord
Throttle_flight = interp1(Thrust_N*n_motors, Throttle, Drag, "pchip", "extrap");
velocity_data = table(V_flight, Range, time_flight, Drag, AspectRatio, Span, ...
    Chord, MassBattery, Power_while_flying, Throttle_flight, CDt, CD_zero, K);

% ---------------- 3. PLOT BLOCK 1 – Range & Endurance --------------------

% ---------------- Figure 1 : Mission‑level performance -------------------
figure;

% --- 1a) Range & Time vs Velocity ----------------------------------------
subplot(1,2,1);                % Place first panel on the left
grid on;                       % Add background grid for readability

yyaxis left                    % LEFT axis → range
plot(V_flight, Range/1000, ...
     'LineWidth', 3, ...
     'Color', [0, 0.4470, 0.7410]);     % Blue curve = mission range
ylabel('Range (km)', ...
       'FontSize', 18, ...
       'FontWeight', 'bold');

yyaxis right                   % RIGHT axis → flight time
plot(V_flight, time_flight, ...
     'LineWidth', 3, ...
     'Color', [0.8500, 0.3250, 0.0980]);% Orange curve = endurance
ylabel('Time of flight (s)', ...
       'FontSize', 18, ...
       'FontWeight', 'bold');

legend('Range (km)', 'Time of flight (s)', ...
       'FontSize', 20, ...
       'FontWeight','normal', ...
       'Location','best');

% Axes cosmetics
ax = gca;
ax.LineWidth = 1.5;
ax.FontSize  = 18;
ax.YAxis(1).Color = 'k';       % Sync left/right axis colours with grid
ax.YAxis(2).Color = 'k';

title('Flight Performance vs Velocity');
xlabel('Flight Speed (m/s)');
% NOTE: zlabel is leftover from an earlier 3‑D version; can be removed if 2‑D.

% --- 1b) Range vs Battery mass & Aspect ratio (contour) ------------------
subplot(1,2,2);                % Second panel on the right

subplot(1,2,2);
battery_mass_mat = repmat(battery_mass, size(D, 1), 1);
contourf(battery_mass_mat, StoreAR, StoreRange/1000, 'LineColor', 'k');
ax = gca; % Get current axes
ax.LineWidth = 1.5; % Make the outer lines bold
ax.FontSize = 18;
grid on;
title('Flight Distance vs Battery Mass', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('Battery Mass (kg)', 'FontSize', 18,'FontWeight', 'bold');
ylabel('Aspect Ratio', 'FontSize', 18, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Range (km)';
cb.Label.FontSize = 20;
cb.FontWeight = 'bold';
grid on;
hold on;
% Overlay contour lines for clarity
[C, h] = contour(battery_mass_mat, StoreAR, StoreRange/1000, 'LineColor', 'black');  % Contour lines
clabel(C, h ,'Color','k','FontSize', 18);  % Label contour lines
grid on;


% ---------------- Figure 2 : Drag vs Velocity vs Aspect Ratio ------------
figure;
plot3(V_flight, Drag, AspectRatio, ...
      'LineWidth', 3, ...
      'Color', [0.8500, 0.3250, 0.0980]);  % Orange curve
grid on;

title('Aspect Ratio Effects on Drag');
xlabel('Flight Speed (m/s)');
ylabel('Drag (N)');
zlabel('Aspect Ratio');

ax = gca;
ax.LineWidth = 1.5;
ax.FontSize  = 18;

end