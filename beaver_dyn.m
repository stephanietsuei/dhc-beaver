function [deriv, ybvel, yuvw, ydl, yfp, ypow, yacc, Caero, Cprop, ...
    FMaero, FMprop, Fgrav, Fwind, yatm, yad1, yad2, yad3] = ...
    beaver_dyn(state, uaero, uprop, uwind)


% Matrices describing the shape of the Beaver aircraft.
AMb = [-0.03554, -0.002226, -0.05504, 0.000591, 0.09448,-0.003117;...
      0.00292,0,-5.578,0,-0.6028,0;...
      5.459,0,0,0,-2.14000000000000,0;...
      -5.16200000000000,0,3.44200000000000,0,0,0;...
      0,-0.767800000000000,0,-0.0618000000000000,0,0.00671900000000000;...
      0,0,0,0,0.692100000000000,0;...
      0,0,0,0,0,0.137300000000000;...
      0,-0.124000000000000,0,-0.504500000000000,0,-0.158500000000000;...
      -0.67480,0,-2.98800000000000,0,-15.5600000000000,0.159500000000000;...
      0,0.3666,0,0.169500000000000,-0.311800000000000,-0.111200000000000;...
      0,0,-0.398000000000000,0,-1.92100000000000,0;...
      -0.0944700000000000,0,-1.37700000000000,0,0.407200000000000,0;...
      0,-0.0295600000000000,0,-0.0991700000000000,0,-0.00387200000000000;...
      0.03412,0.115800000000000,0,0.00693400000000000,0,-0.0826500000000000;...
      1.10600000000000,0,-1.26100000000000,0,0,0;...
      0,0.523800000000000,0,0,0,0;...
      0,0,0,-0.0826900000000000,0,0;...
      0,0,-15.9300000000000,0,0,0;...
      0,-0.160000000000000,0,0,0,0;];
AM = AMb';
EM = [0.116100000000000,0,0.145300000000000,0;...
      0,0,0,0;...
      -0.156300000000000,0,0,0;...
      0,0,0,-0.0140600000000000;...
      -0.0789500000000000,0,0,0;...
      0,-0.00302600000000000,0,0;];
GM1 = [1.58750000000000,14.6300000000000,23.2300000000000,5368.39000000000,...
    6928.93000000000,11158.7500000000,0,117.640000000000,0,2288.23100000000;];
GM2 = [0.000186318630694235,0,1.96424543204838e-06,...
                0,0.0188532401483411,0,0,-0.788325344315717,0;...
       0,0.000144322427849610,0,-0.0169780904122282,...
                0,0.835678813323269,0,0,0.0169780904122282;...
       1.96424543204838e-06,0,8.96364802359248e-05,0,...
                -0.139650239034744,0,0,-0.0188532401483411,0;];

            
% Variables defined as part of GM1, GM2, AM, and EM
cbar = GM1(1);
b = GM1(2);
S = GM1(3);
m = GM1(10);
grav = 9.80665;
CYbdot = AM(2,19);


% number of timesteps we're computing derivatives for at once
num_timesteps = size(state,2);


% Compute the derivative and the output vector
[yatm, yad1, yad2, yad3] = airdata(state);
[ydl, Caero, FMaero] = aerodynamics(state, uaero, yad1);
[ypow, Cprop, FMprop] = engine(state, uprop, yatm, yad1);
Fgrav = gravity(state, yatm);
Fwind = fwind(state, uwind);
[Ftot, Mtot] = fmsort(FMaero, FMprop, Fgrav, Fwind);
yhlp = hlpfcn(state);
[deriv, ybvel] = eqns(state, Ftot, Mtot, uwind, yatm, yhlp);
[yfp, yuvw, yacc] = addlouts(state, deriv, yhlp, Ftot, Fgrav);



%------------------------------------------------------------------------------%
% All the nested functions! Copied directly from the original Simulink model.  %
%------------------------------------------------------------------------------%


% Computes the values of atmospheric variables that vary with the state of
% the plane.
function [yatm, yad1, yad2, yad3] = airdata(x)

    yatm = atmosph(x);
    yad1 = airdata1(x, yatm);
    yad2 = airdata2(yatm, yad1);
    yad3 = airdata3(x, yatm, yad1);

    % Computes the air density, pressure, temperature, gravity
    function yatm = atmosph(x)
        g = grav*(6371020 ./ (6371020 + x(12,:))).^2;
        T = 288.15 - .0065*x(12,:);

        u1 = [g; T];

        ps = 101325 * (u1(2,:) / 288.15).^(u1(1,:) / 1.86584);
        mu = (1.458*10^(-6) * u1(2,:).^1.5) ./ (u1(2,:) + 110.4);

        u2 = [u1; ps];

        rho = u2(3,:) ./ (287.053 * u2(2,:));

        yatm = [rho; ps; T; mu; g];
    end

    % Computes the speed of sound (a), the Mach number (M), and the dynamic
    % pressure (qdyn).
    function yad1 = airdata1(x, yatm)
        a = sqrt(401.8743 * yatm(3,:));

        u_M = [a; x];
        M = u_M(2,:) ./ u_M(1,:);

        u_q = [x; yatm];
        qdyn = .5 * u_q(13,:) .* u_q(1,:).^2;

        yad1 = [a; M; qdyn];
    end

    % Computes the impact pressure (qc), the equivalent airspeed (Ve), and
    % the calibrated airspeed (Vc)
    function yad2 = airdata2(yatm, yad1)
        u = [yatm; yad1];

        qc = u(2,:) .* ((1 + .2*u(7,:).^2).^3.5 - 1);
        Ve = sqrt(u(8,:) * 1.63265);

        Vc = sqrt(579000*((1 + qc(1,:)/101325).^(1/3.5)-1));

        yad2 = [qc; Ve; Vc];
    end

    % Computes the total Temperature (Tt), and the Reynolds Numbers (Re,
    % Rc)
    function yad3 = airdata3(x, yatm, yad1)

         % loaded at beginning of file

        u = [x; yatm; yad1];

        Tt = u(15,:) .* (1 + .2*u(19,:).^2);
        Re = (u(13,:).*u(1,:)) ./ u(16,:);
        Rc = cbar*Re;

        yad3 = [Tt; Re; Rc];
    end
end


% Computes aerodynamic forces and moments for the aircraft
function [ydl, Caero, FMaero] = aerodynamics(x, uaero, yad1)

    ydl = dimless(x);
    Caero = aeromod(x, uaero, ydl);
    FMaero = fmdims_aero(yad1, Caero);

    % Computes the rotational velocities. This is the reason we cannot have
    % zero velocity.
    function ydl = dimless(x)
        pb2V = .5 * x(4,:) * b ./ x(1,:);
        qcV = x(5,:) * cbar ./ x(1,:);
        rb2V = .5 * x(6,:) * b ./ x(1,:);

        ydl = [pb2V; qcV; rb2V];
    end

    % Defines aerodynamic force and moment coefficients
    function Caero = aeromod(x, uaero, ydl)
        alpha = x(2,:);
        beta = x(3,:);
        deltae = uaero(1,:);
        deltaa = uaero(2,:);
        deltar = uaero(3,:);
        deltaf = uaero(4,:);
        
        u1 = zeros(7,num_timesteps);
        u1(1,:) = ones(1,num_timesteps);
        u1(2,:) = alpha;
        u1(3,:) = alpha.^2;
        u1(4,:) = alpha.^3;
        u1(5,:) = beta;
        u1(6,:) = beta.^2;
        u1(7,:) = beta.^3;

        u2 = zeros(9,num_timesteps);
        u2(1,:) = deltae;
        u2(2,:) = deltaf;
        u2(3,:) = deltaa;
        u2(4,:) = deltar;
        u2(5,:) = alpha.*deltaf;
        u2(6,:) = alpha.*deltar;
        u2(7,:) = alpha.*deltaa;
        u2(8,:) = beta.^2 .* deltae;
        u2(9,:) = zeros(1,num_timesteps);

        u = [u1; ydl; u2];

        Caero = AM * u;
    end

    % Computes forces and moments on the airplane
    function FMaero = fmdims_aero(yad1, Caero)
        u = [yad1; S*Caero];

        FMaero = zeros(6,num_timesteps);
        FMaero(1,:) = u(4,:).*u(3,:);
        FMaero(2,:) = u(5,:).*u(3,:);
        FMaero(3,:) = u(6,:).*u(3,:);
        FMaero(4,:) = b*u(7,:).*u(3,:);
        FMaero(5,:) = cbar*u(8,:).*u(3,:);
        FMaero(6,:) = b*u(9,:).*u(3,:);
    end

end


% Computes forces and coefficients from the plane's engine
function [ypow, Cprop, FMprop] = engine(x, uprop, yatm, yad1)

    ypow = power(x, uprop, yatm);
    Cprop = engmod(x, ypow);
    FMprop = fmdims_prop(yad1, Cprop);
    
    function ypow = power(x, uprop, yatm)
        u1 = [x; uprop; yatm];
        P = .7355*(-326.5 + .00412*(u1(14,:)+7.4).*(u1(13,:)+2010) + ...
            (408-.0965*u1(13,:)).*(1-u1(15,:)/1.225));
        
        u2 = [u1; P];
        dpt = .08696 + 191.18*(u2(20,:)./(.5*u2(15,:).*u2(1,:).^3));
        
        ypow = [dpt; P];
    end


    function Cprop = engmod(x, ypow)
        dpt = ypow(1,:);
        alpha = x(2,:);
        
        u = [dpt; alpha];
        ytmp = [dpt; ...
                u(1,:).^3; ...
                u(2,:).*u(1,:).^2;
                u(1,:).*u(2,:).^2;];
        Cprop = EM*ytmp;
    end


    function FMprop = fmdims_prop(yad1, Cprop)
        u = [yad1; S*Cprop];

        FMprop = zeros(6,num_timesteps);
        FMprop(1,:) = u(4,:).*u(3,:);
        FMprop(2,:) = u(5,:).*u(3,:);
        FMprop(3,:) = u(6,:).*u(3,:);
        FMprop(4,:) = b*u(7,:).*u(3,:);
        FMprop(5,:) = cbar*u(8,:).*u(3,:);
        FMprop(6,:) = b*u(9,:).*u(3,:);
    end
end


function Fgrav = gravity(x, yatm)
    u = [x; yatm];
    
    Xgr = -m*u(17,:).*sin(u(8,:));
    Ygr = m*u(17,:).*cos(u(8,:)).*sin(u(9,:));
    Zgr = m*u(17,:).*cos(u(8,:)).*cos(u(9,:));
    
    Fgrav = [Xgr; Ygr; Zgr];
end


% Computes forces from wind disturbances
function Fwind = fwind(x, uwind)
    u = [x; uwind];
    
    Xwm = u(16,:) + u(5,:).*u(15,:) - u(6,:).*u(14,:);
    Ywm = u(17,:) - u(4,:).*u(15,:) + u(6,:).*u(13,:);
    Zwm = u(18,:) + u(4,:).*u(14,:) - u(5,:).*u(13,:);
    
    Fwind = -m*[Xwm; Ywm; Zwm];
end


% Sums up all the different forces and moments
function [Ftot, Mtot] = fmsort(FMaero, FMprop, Fgrav, Fwind)
    Ftot = FMaero(1:3,:) + FMprop(1:3,:) + Fgrav + Fwind;
    Mtot = FMaero(4:6,:) + FMprop(4:6,:);
end


function [xdot, ybvel] = eqns(x, Ftot, Mtot, uwind, yatm, yhlp)
    ybvel = uvw(yhlp, x);
    ybvel2 = ybvel + uwind(1:3,:);
    xdot = odes(x, Ftot, Mtot, yhlp, ybvel2);
    xdot = xdotcorr(xdot, yhlp, yatm);
    
    function ybvel = uvw(yhlp, x)
        u = [x; yhlp];
        u1 = u(1,:).*u(13,:).*u(15,:);
        v = u(1,:).*u(16,:);
        w = u(1,:).*u(14,:).*u(15,:);
        
        ybvel = [u1; v; w];
    end
    
    
    % Computes the state derivative.
    function xdot = odes(x, Ftot, Mtot, yhlp, ybvel2)
        u = [x; Ftot; Mtot; yhlp];
        dot1 = Vabdot(u);
        dot2 = pqrdot(u);
        dot3 = eulerdot(u);
        dot4 = xyhdot(ybvel2, u);
        xdot = [dot1; dot2; dot3; dot4];
        
        % Time derivative of the speed, angle of attack, and angle of
        % sideslip
        function dot1 = Vabdot(u)
            Vdot = (u(13,:).*u(19,:).*u(21,:) + u(14,:).*u(22,:) + ...
                u(15,:).*u(20,:).*u(21,:)) / m;
            adot = (-u(13,:).*u(20,:) + u(15,:).*u(19,:)) ./ ...
                   (m*u(1,:).*u(21,:)) - ...
                   u(23,:).*(u(4,:).*u(19,:) + u(6,:).*u(20,:)) + u(5,:);
            bdot = (-u(13,:).*u(19,:).*u(22,:) + u(14,:).*u(21,:) ...
                - u(15,:).*u(20,:).*u(22,:)) ./ ...
                (m*u(1,:)) + u(4,:).*u(20,:) - u(6,:).*u(19,:);
            dot1 = [Vdot; adot; bdot];
        end
        
        % Time derivative of rotational velocities
        function dot2 = pqrdot(u)
            pqr = u(4:6,:);
            lmn = u(16:18,:);
            
            ytmp = [lmn; ...
                    pqr(1,:).^2; ...
                    pqr(1,:).*pqr(2,:); ...
                    pqr(1,:).*pqr(3,:); ...
                    pqr(2,:).^2; ...
                    pqr(2,:).*pqr(3,:); ...
                    pqr(3,:).^2];
            dot2 = GM2*ytmp;
        end
        
        % Time derivative of attitude angles
        function dot3 = eulerdot(u)
            psidot = (u(5,:).*u(28,:) + u(6,:).*u(29,:)) ./ u(27,:);
            thetadot = u(5,:).*u(29,:) - u(6,:).*u(28,:);
            up = [u; psidot];
            phidot = up(4,:) + up(30,:).*up(26,:);
            
            dot3 = [psidot; thetadot; phidot];
        end
        
        % Time derivative of position coordinates
        function dot4 = xyhdot(ybvel2, u)
            u1 = [u; ybvel2];
            
            out11 = u1(30,:).*u1(27,:) + (u1(31,:).*u1(28,:) + ...
                u1(32,:).*u1(29,:)).*u1(26,:);
            out12 = u1(31,:).*u1(29,:) - u1(32,:).*u1(28,:);
            
            u2 = [u; out11; out12];
            xedot = u2(30,:).*u2(25,:) - u2(31,:).*u2(24,:);
            yedot = u2(30,:).*u2(24,:) + u2(31,:).*u2(25,:);
            hdot = u1(30,:).*u1(26,:) - (u1(31,:).*u1(28,:) + ...
                u1(32,:).*u1(29,:)).*u1(27,:);
            
            dot4 = [xedot; yedot; hdot];
        end
    end
    
    function xdcorr = xdotcorr(xdot, yhlp, yatm)
        u = [yhlp; xdot(3,:); yatm];
        xdcorr = xdot;
        xdcorr(3,:) = u(12,:) ./ (1 - (u(3,:).*u(13,:)*b*S*CYbdot)/(4*m));
    end
end



function [yfp, yuvw, yacc] = addlouts(x, xdot, yhlp, Ftot, Fgrav)

    yfp = flpath(x, xdot, yhlp);
    yuvw = uvwdot(x, xdot, yhlp);
    yacc = accel(Ftot, Fgrav);
    
    % The flight path angle
    function yfp = flpath(x, xdot, yhlp)
        u = [x; xdot; yhlp];
        gamma = asin(u(24,:)./u(1,:));
        fpa = u(13,:) / grav;
        chi = u(3,:) + u(7,:);
        Phi = asin(u(34,:).*u(33,:));
        
        yfp = [gamma; fpa; chi; Phi];
    end
    
    function yuvw = uvwdot(x, xdot, yhlp)
        u = [yhlp(1:4,:); x(1,:); xdot];
        
        udot = u(6,:).*u(1,:).*u(3,:) - u(5,:).*(u(7,:).*u(2,:).*u(3,:) + ...
            u(8,:).*u(1,:).*u(4,:));
        vdot = u(6,:).*u(4,:) + u(5,:).*u(3,:).*u(8,:);
        wdot = u(6,:).*u(2,:).*u(3,:) + u(5,:).*(u(7,:).*u(1,:).*u(3,:) ...
            - u(8,:).*u(2,:).*u(4,:));
        
        yuvw = [udot; vdot; wdot];
    end
    
    function yacc = accel(Ftot, Fgrav)
        u = [Ftot; Fgrav];
        weight = m*grav;
        
        Ax = (u(1,:)-u(4,:)) / weight;
        Ay = (u(2,:)-u(5,:)) / weight;
        Az = (u(3,:)-u(6,:)) / weight;
        axk = u(1,:) / weight;
        ayk = u(2,:) / weight;
        azk = u(3,:) / weight;
        
        yacc = [Ax; Ay; Az; axk; ayk; azk];
    end
end


% Computes sines and cosines to make computing the equations of motion
% easier
function yhlp = hlpfcn(x)
    yhlp = [cos(x(2,:)); ...
            sin(x(2,:)); ...
            cos(x(3,:)); ...
            sin(x(3,:)); ...
            tan(x(3,:)); ...
            sin(x(7,:)); ...
            cos(x(7,:)); ...
            sin(x(8,:)); ...
            cos(x(8,:)); ...
            sin(x(9,:)); ...
            cos(x(9,:))];
end



end