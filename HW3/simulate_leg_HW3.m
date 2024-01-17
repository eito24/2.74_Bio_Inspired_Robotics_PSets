function simulate_leg_HW3()

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
    
    %% Parameter vector
    p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 10;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [-pi/4; pi/2; 0; 0];
    z_out = zeros(4,num_step);
    z_out(:,1) = z0;
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p);
        z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;
        z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt;
    end

    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    rE = zeros(3,length(tspan));
    for i = 1:length(tspan)
        rE(:,i) = position_foot(z_out(:,i),p);
    end
    figure(2); clf;
    plot(tspan,rE(1:2,:))
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','y'});

    
    %% Animate Solution
    figure(3); clf;
    hold on
   
    %% Optional, plot foot target information
    
    % target position Q1.4
    plot(.025 , -.125, 'r.','MarkerSize',6); 
    
    % Target traj. Q 1.6
    plot( .025*cos(0:.01:2*pi), -.125+.025*sin(0:.01:2*pi),'k--'); 
    
    % Ground Q2.3
    plot([-.2 .2],[-.125 -.125],'k'); 
    
    animateSol(tspan,z_out,p);
end

function tau = control_law(t,z,p)
    % t is current time, z=[th1,th2,dth1,dth2], p is parameters
    % Controller gains, Update as necessary for Problem 1
    K_x = 40; % Spring stiffness X
    K_y = 40; % Spring stiffness Y
    D_x = 4;  % Damping X
    D_y = 4;  % Damping Y
    w = 3;
    drEd = [-0.025*w*sin(w*t) 0.025*w*cos(w*t)];
    rEd = [0.025*cos(w*t) -0.125+0.025*sin(w*t)]'; % Desired position of foot
    
    %% STEPS TO COMPLETE PROBLEM 1.3
    % a. Compute r_E
    % b. Compute J, the jacobian of r_E
    % c. Use these results to compute \tau as specified in the write-up
    % J = dr/dq
    % dV/q = dV/dr * dr/dq = dV/dr * J
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.
    J = jacobian_rE(z,p);
    J = J(1:2,1:2);
    x_d = rEd(1);
    y_d = rEd(2);
    rE = position_foot(z,p);
    x = rE(1);
    y = rE(2);
    dr_E = drE(z,p);
    dx = dr_E(1);
    dy = dr_E(2);
    F = [K_x*(x-x_d); K_y*(y-y_d)]+[D_x*(dx-drEd(1)); D_y*(dy-drEd(2))];
    tau = -J'*F;
end


function Fc = contact_force(z,p)

    %% Fixed parameters for contact
    K_c = 100;
    D_c = 2;
    yC  = -.125;
    rE = position_foot(z,p);
    dr_E = drE(z,p);
    dy = dr_E(2);
    y = rE(2);
    %% STEPS TO COMPLETE PROBLEM 2.1
    % a. Compute constraint C which gives height of foot relative to ground
    % b. Compute constraint rate, \dot{C}
    % c. Set Fc based on compliant contact model
    % d. If foot is above the ground, or Fc<0, set Fc = 0
    C = y-yC;
    dC = dy;
    % Hint: Some of the functions generated by derive_leg will help with 
    % steps a and b.
    if C > 0
        Fc = [0;0;];
    else
        Fc = [0; -K_c*C - D_c*dC];
    end
 % replace this line using steps a-d to compute Fc
end


function Tauc = joint_limit_torque(z,p)
    %% Fixed parameters for rotational spring damper at joint limit contact
    Kappa_c = 10;
    Dampa_c = 0.2;

    Tauc = 0;
end


function dz = dynamics(t,z,p)
    th1 = z(1);     th2 = z(2);
    dth1= z(3);     dth2= z(4);
    
    % Get mass matrix
    A = A_leg(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_leg(z,tau,p);
  
    % Compute the contact force (used for problem 2)
    Fc = contact_force(z,p);
   
    % Compute the contribution of the contact force to the generalied force
    J = jacobian_rE(z,p);
    J = J(1:2,1:2);
    QFc = J'*Fc;  %% YOUR CODE HERE for Q2.2

    % Compute the contact force (used for problem 2.5)
    Tauc = joint_limit_torque(z,p);
    QTauc= [0; 0];
    
    % Solve for qdd.
    qdd = A\(b + QFc + QTauc);
    dz = 0*z;
    
    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
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