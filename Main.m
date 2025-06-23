% =========================================================================
%  HERMES eVTOL / Tailsitter Concept – Performance Optimisation Script
% =========================================================================
%  The aerodynamics have been pre-corrected for a **biplane** configuration.
%  ------------------------------------------------------------------------
%  REQUIRED DATA FILES
%    • airfoil_data.mat — provides struct `airfoil_e422` with fields:
%        - alpha : angle of attack vector (deg)
%        - CL    : lift-coefficient vector (dimensionless)
%    • MADMotor.mat    — provides struct `MADMotor` with fields:
%        - Power_W : electrical power vs. static thrust (W)
%        - Thrust_N: static thrust values (N)
%        - RPM     : propeller speed (rev · min⁻¹)
% (If the airfoil or motor changed make sure to load it in the save way)
%  ------------------------------------------------------------------------
%  OUTPUT VARIABLES
%    • CL_data, bestCL         – optimum lift-coefficients
%    • velocity_data           – optimum cruise velocity at each design point
%    • diff_velocity_data      – margin between optimum and limit velocity
%  ------------------------------------------------------------------------
%  DEPENDENCIES (user-supplied m-files)
%    compute_vmax_prop.m : Prop-tip Mach / J-limit speed calculator
%    OptFun.m            : CL optimiser
%    VelOpt.m            : Velocity optimiser
%  ------------------------------------------------------------------------
%  Author : Bartolome Cabrer Falomir
%  Created: 17-Jun-2025 (Europe/London)
%  License: Imperial College London undergraduate license
% =========================================================================

% --------------------------------------------------------------
% 📌 CONFIG NOTE – Switching from X-wing/biplane to single wing
%
%  ▸ No X-wing geometry?
%       thetaDeg  = 0;        % remove the dihedral break
%       length_x  = 0;        % no in-board fold section
%
%  ▸ Only one wing plane (not a biplane)?
%       b = b / 2;            % halve the span value
%  
%       Cross_section_N = Cross_section_N/2; % halve the number of spars
%
%       (Because the code accounts for 2 wings so if not will doblue the 
%        number of ribs) 
%
%  Doing the three edits above turns the model into a simple,
%  single-deck, flat wing without the X-tail “kink”.
% --------------------------------------------------------------
%% 1) House-keeping --------------------------------------------------------
% clear;

%% 2) Load aerodynamic & propulsion data ----------------------------------
load("airfoil_data.mat")   % ⇢ loads struct `airfoil_e422`
load("MADMotor.mat")       % ⇢ loads struct `MADMotor`

% Convenience aliases
Power_W_raw  = MADMotor.Power_W;  % Electrical power draw (W)
Thrust_N_raw = MADMotor.Thrust_N;  % Corresponding static thrust (N)
RPM_raw      = MADMotor.RPM;      % Propeller speed (rev/min)

%% 3) Build a uniform throttle grid ---------------------------------------
a_thrust   = linspace(0, max(Thrust_N_raw), 20);   % [N] 20-point thrust vector
Throttle   = a_thrust / max(a_thrust);             % [-] Normalised throttle 0-1

% Interpolate electrical power (monotone PCHIP) and over-estimate at low thrust
Power_grid = interp1(Thrust_N_raw, Power_W_raw, a_thrust, "pchip", "extrap");
Power_grid(Power_grid < min(Power_W_raw)) = min(Power_W_raw);

% Interpolate RPM on the same grid
RPM_grid   = interp1(Thrust_N_raw, RPM_raw, a_thrust, "pchip", "extrap");

% Column vectors to match downstream functions
Power_W  = Power_grid.';      % [W]
Thrust_N = a_thrust.';        % [N]
RPM      = RPM_grid.';        % [rev/min]

%% 4) Propeller geometry & speed limits -----------------------------------
a_sound     = 343;       % [m/s] speed of sound at 20 °C
D           = ;     % [m] propeller diameter ( **CHANGE HERE** if needed )
J_max       = 0.8;       % [-] max advance ratio before stall
M_tip_max   = 0.75;      % [-] max allowable tip Mach number

max_speed = compute_vmax_prop(D, max(RPM), J_max, M_tip_max, a_sound);

% Optional: power available using left-continuous interpolation (was `pa` in
% the original script).  Not used in OptFun/VelOpt but kept for reference.
P_available = interp1(Throttle, Power_W, Throttle, "previous", "extrap");

%% 5) Fleet configuration & physical constants ----------------------------
n_motors   = ;          % Number of propulsive units
rho        = 1.2;        % [kg/m³] ISA sea-level air density
mu         = 1.82e-5;    % [Pa·s]  air dynamic viscosity (for Re calc.)

% Operational flight envelope
max_safety_V   = ;     % [m/s] hard safety velocity limit
min_velocity   = ;     % [m/s] lowest considered flight speed
velocity_step  = 1;      % [m/s] resolution
V_flight     = min_velocity:velocity_step:max_safety_V;  % speed vector

% Thrust limit (at MT % throttle) for envelope calculations
MT_percent  = 80;  % [%] throttle cap for steady-state cruise sizing
Thrust_max  = interp1(Throttle, Thrust_N .* n_motors, MT_percent, "pchip", "extrap");

%% 6) Mass breakdown -------------------------------------------------------
Fuselage_mass   = ;  % [kg] structure + avionics + etc.
Motor_mass    = ;            % [kg] per motor
Payload_mass  = ;           % [kg] nominal human payload
Dry_mass      = Fuselage_mass + Payload_mass + Motor_mass * n_motors;

% Battery mass sweep and corresponding capacities
battery_mass      = 50:5:80;                 % [kg] design sweep
battery_specificE = ;   % [J/kg] 230 Wh/kg in SI
battery_specificE =  battery_specificE * 3600;             
battery_capacity  = battery_specificE .* battery_mass;  % [J]

% Total aircraft weight (incl. battery) in Newtons
Weight_tot = (Dry_mass + battery_mass) * 9.81;   % [N]

%% 7) Airfoil data pre‑processing ----------------------------------------
CL_raw  = airfoil_e422.CL;     % dimensionless lift coefficient vector
AoA_raw = airfoil_e422.alpha;  % [deg] matching angles of attack

% Use only CL ≤ 0.45 for reliable linear region (indices 21:73)
CL  = CL_raw(21:73);
AoA = AoA_raw(21:73);
%% 8) Wing & fuselage geometry --------------------------------------------
% Fuselage
fuse_diam = ;   % [m] diameter
fuse_len  = ;     % [m] cylindrical length considered in wetted area

length_x = ; % [m] horiz. projection of X. From fuselage to wing
max_span = ; %[m] the max available span. Tip to tip without fuselage

% Wing planform
b_span   = (2 * length_x : 0.1 : max_span)';  % [m] wingspan design vector
c_min    = ;                 % [m] minimum allowable chord

% Aerodynamic correction factors
thetaDeg = ;       % [deg] X wing angle
kM       = 1;      % [-] span efficiency factor (≈ 1 for good twist)
Horner_f = 0.0119; % [-] Horner factor (1951) for zero sweep (Cd₀ corr.)



%% 9) Airfoil geometry for structural sizing ------------------------------
% Interpolate upper & lower surface coordinates to a common x grid for
% later structural weight estimation.  Simply assume `top` and `bottom`
% variables already exist in workspace (otherwise load them here).

dx  = 0.005;        % step in x (0–1) for interpolation
xq  = 0:dx:1;       % common abscissa

% Ensure top / bottom spline have been defined before this block!
top_interp    = [xq', interp1(top(:,1),    top(:,2),    xq)'];
bottom_interp = [xq', interp1(bottom(:,1), bottom(:,2), xq)'];

% Structural discretisation & material props
geom.top              = top_interp;   % upper surface coords
geom.bottom           = bottom_interp; % lower surface coords
geom.length_x         = length_x;      % [m] horiz. projection of dihedral
geom.xq               = xq;           % shared chordwise stations
geom.dx               = dx;           % spacing
geom.Cross_section_N  = ;            % # ribs / wing segment
geom.Spar_length      = ;         % [m] spar square side
geom.aero_thickness   = ;         % [m] profile thickness for rib 
geom.skin_thickness   = ;       % [m] composite skin
geom.rho_spar         = ;         % [kg/m³] density of spars material
geom.rho_skin         = ;         % [kg/m³] CFRP skin
geom.rho_rib          = ;         % [kg/m³] density of ribs material
geom.thetaDeg         = thetaDeg;     % pass dihedral to functions

%% 10) Optimise lift coefficient schedule ---------------------------------
[CL_data, bestCL] = OptFun(Power_W, Thrust_N, Throttle, n_motors, ...
    CL, V_flight, Weight_tot, battery_mass, b_span, battery_capacity, ...
    rho, c_min, fuse_len, Horner_f, kM, fuse_diam, mu, Thrust_max, AoA, ...
    thetaDeg, geom, max_speed, geom.length_x, Payload_mass);

%% 11) Compute optimum velocity envelope ----------------------------------
[velocity_data, diff_velocity_data] = VelOpt(Power_W, Thrust_N, Throttle, n_motors, ...
    bestCL, V_flight, Weight_tot, battery_mass, b_span, battery_capacity, ...
    rho, c_min, fuse_len, Horner_f, kM, fuse_diam, mu, Thrust_max, RPM, ...
    thetaDeg, geom, max_speed, geom.length_x);

%% 12) (Optional) Save processed data -------------------------------------
% Uncomment the next line to save variables for post-processing / plotting
% save('Hermes_data.mat','velocity_data','CL_data','diff_velocity_data')

%% End of script ----------------------------------------------------------
