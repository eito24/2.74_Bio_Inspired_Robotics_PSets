function simulate_leg()
    %% Definte fixed paramters
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.81;    
    
    restitution_coeff = 0.;
    friction_coeff = 10;
    ground_height = -0.125;
    %% Parameter vector
    p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';
       
    %% Simulation Parameters Set 2 -- Operational Space Control
    p_traj.omega = 3;
    p_traj.x_0   = 0;
    p_traj.y_0   = -.125;
    p_traj.r     = 0.025;
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 10;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [-pi/4; pi/2; 0; 0];
    z_out = zeros(4,num_step);
    z_out(:,1) = z0;
    
    for i=1:num_step-1
        z = z_out(:,i);
        dz = dynamics(tspan(i), z, p, p_traj);
        % Velocity update with dynamics
        new_z = z_out(:,i) + dz*dt;
        qdot = discrete_impact_contact(new_z, p, restitution_coeff, friction_coeff, ground_height);
        z_out(:,i+1) = [new_z(1:2);qdot];
        % Position update
        z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt;
    end
    
    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    rE = zeros(2,length(tspan));
    vE = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_foot(z_out(:,i),p);
        vE(:,i) = velocity_foot(z_out(:,i),p);
    end
    
    figure(2); clf;
    plot(tspan,rE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,rE(2,:),'b','LineWidth',2)
    plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
    
    
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','x_d','y','y_d'});

    figure(3); clf;
    plot(tspan,vE(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,vE(2,:),'b','LineWidth',2)
    
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x','vel_y'});
    
    figure(4)
    plot(tspan,z_out(1:2,:)*180/pi)
    legend('q1','q2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    
    figure(5)
    plot(tspan,z_out(3:4,:)*180/pi)
    legend('q1dot','q2dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    
    %% Animate Solution
    figure(6); clf;
    hold on
   
  
    % Target traj
    TH = 0:.1:2*pi;
    plot( p_traj.x_0 + p_traj.r * cos(TH), ...
          p_traj.y_0 + p_traj.r * sin(TH),'k--'); 
    
    % Ground Q2.3
    plot([-.2 .2],[ground_height ground_height],'k'); 
    
    animateSol(tspan, z_out,p);
end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 150.; % Spring stiffness X
    K_y = 150.; % Spring stiffness Y
    D_x = 10.;  % Damping X
    D_y = 10.;  % Damping Y

    % Desired position of foot is a circle
    omega_swing = p_traj.omega;
    rEd = [p_traj.x_0 p_traj.y_0 0]' + ...
            p_traj.r*[cos(omega_swing*t) sin(omega_swing*t) 0]';
    % Compute desired velocity of foot
    vEd = p_traj.r*[-sin(omega_swing*t)*omega_swing    ...
                     cos(omega_swing*t)*omega_swing   0]';
    % Desired acceleration
    aEd = p_traj.r*[-cos(omega_swing*t)*omega_swing^2 ...
                    -sin(omega_swing*t)*omega_swing^2 0]';
    
    % Actual position and velocity 
    rE = position_foot(z,p);
    vE = velocity_foot(z,p);
    
    % Compute virtual force for Question 1.4 and 1.5
    J = jacobian_foot(z,p);
    Jdot = jacobian_dot_foot(z,p);
    M = A_leg(z,p);
    V = Corr_leg(z,p);
    G = Grav_leg(z,p);

    effective_mass = inv(J*inv(M)*J');
    coriolis_centripetal_force = effective_mass*J*inv(M)*V-effective_mass*Jdot*z(3:4);
    gravity_force = effective_mass*J*inv(M)*G;
    f  = [K_x * (rEd(1) - rE(1) ) + D_x * ( - vE(1) ) ;
          K_y * (rEd(2) - rE(2) ) + D_y * ( - vE(2) ) ];
    newf_1_2 = effective_mass*(aEd(1:2,:)+f)+coriolis_centripetal_force+gravity_force;
    newf_no_g = effective_mass*(aEd(1:2,:)+f)+coriolis_centripetal_force;
    newf_no_c = effective_mass*(aEd(1:2,:)+f)+gravity_force;
    newf_no_a = effective_mass*f+coriolis_centripetal_force+gravity_force;
    
    %% Task-space compensation and feed forward for Question 1.8

    % Map to joint torques  
    tau = J' * newf_1_2;
end


function dz = dynamics(t,z,p,p_traj)
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
    
    % Solve for qdd.
    qdd = A\(b); % [ddth1 ddth2]
    dz = 0*z;
    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)
    rE = position_foot(z,p);
    y = rE(2);
    C_y = y-yC;
    vE = velocity_foot(z,p);
    ydot = vE(2);
    Cdoty = ydot;
    J = jacobian_foot(z,p);
    M = A_leg(z,p);
    effective_mass_y = inv(J(2,:)*inv(M)*J(2,:)');
    effective_mass_x = inv(J(1,:)*inv(M)*J(1,:)');
    qdot = z(3:4);
    if C_y<0 && Cdoty<0
        vert_impulse_force = effective_mass_y*(-rest_coeff*Cdoty-J(2,:)*z(3:4));
        qdot = qdot + inv(M)*J(2,:)'*vert_impulse_force;
        tangential_impulse_force = effective_mass_x*(0-J(1,:)*qdot);
        if tangential_impulse_force > fric_coeff*vert_impulse_force
            tangential_impulse_force = fric_coeff*vert_impulse_force;
        elseif tangential_impulse_force < -fric_coeff*vert_impulse_force
            tangential_impulse_force = -fric_coeff*vert_impulse_force;
        end
        qdot = qdot + inv(M)*J(1,:)'*tangential_impulse_force;
    end
end

function qdot = joint_limit_constraint(z,p)

end

function animateSol(tspan, x,p)
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 .2 -.3 .1]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_leg(z,p);

        rA = keypoints(:,1); % Vector to base of cart
        rB = keypoints(:,2);
        rC = keypoints(:,3); % Vector to tip of pendulum
        rD = keypoints(:,4);
        rE = keypoints(:,5);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[0 rB(1)]);
        set(h_OB,'YData',[0 rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);

        pause(.01)
    end
end