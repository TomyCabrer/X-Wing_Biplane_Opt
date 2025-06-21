function [new_c,F] = chordTesting(c, b, ConstWeight, Vmin_lift, CL_loop, rho, geom)
    % Calculate current wing mass for this c and b
    %WHEN THIS WORKS USE PREVIOUS RESULT AS GUESS
    w_wing = wing_weight_finder( ...
                 geom.top, geom.bottom, b, c, ...
                 geom.min_b, geom.xq, geom.dx)

    % Total weight
    Wtot = ConstWeight + w_wing;

    % Required chord from the lift equation
    c_req = Wtot / (rho * Vmin_lift^2 * CL_loop * b)

    % Residual (we want this = 0)
    
    F =  c_req-c;
    new_c=c-F;
end


        V_flight(n, 1) = V(1);
        time_flight(n, 1) = tf(idx) / 60; % Convert flight time to minutes
        Range(n, 1) = R(idx);
        Drag(n, 1) = D(idx);
        AspectRatio(n, 1) = AR(idx);
        MassBattery(n, 1) = battery_mass(col);
        Power_while_flying(n, 1) = power_flight(idx);
        CDt(n, 1) = CD(idx);
        CD_zero(n ,1) = CD0(idx);
        K(n ,1) = k(idx);

                V_flight(m, n) = V(1);
        time_flight(m, n) = tf(idx) / 60; % Convert flight time to minutes
        Range(m, n) = R(idx);
        Drag(m, n) = D(idx);
        AspectRatio(m, n) = AR(idx);
        Span(m, n) = b(row);
        MassBattery(m, n) = battery_mass(col);
        Power_while_flying(m, n) = power_flight(idx);
        CDt(m, n) = CD(idx);
        VL(m, n) = Vmin_lift;