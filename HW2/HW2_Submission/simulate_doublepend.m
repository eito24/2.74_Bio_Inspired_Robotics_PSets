function simulate_doublepend()

    %% Definte fixed paramters
    l1 = 1;
    l2 = 0.5;
    c1 = 0.5;
    c2 = 0.25;
    m1 = 1;
    m2 = 1;
    I1 = 0.05;
    I2 = 0.05;
    g = 9.81;

    p   = [c1; c2; l1; l2; m1; m2; I1; I2; g;];        % parameters

    %% Perform Dynamic simulation    
    dt = 0.001;
    tf = 10;
    num_steps = floor(tf/dt);
    tspan = linspace(0, tf, num_steps); 
    z0 = [3; 0; 0; 0];
    z_out = zeros(4,num_steps);
    z_out(:,1) = z0;
    for i=1:num_steps-1
        dz = dynamics(z_out(:,i), p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
    end
    final_state = z_out(:,end);
    

    %% Compute Energy
    E = energy_doublepend(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');

    %% Animate Solution
    figure(2); clf;
        % Prepare plot handles
    hold on
    h_pole1 = plot([0],[0],'LineWidth',3);
    h_pole2 = plot([0],[0],'LineWidth',3);
    xlabel('x')
    ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-2 2 -2 1]);
    skip_frame = 10;

    %Step through and update animation
    for i=1:num_steps
        if mod(i, skip_frame)
            continue
        end
        % interpolate to get state at current time.
        t = tspan(i);
        z = z_out(:,i);
        keypoints = keypoints_doublepend(z,p);

        rB = keypoints(:,2); % Vector to tip of rod 1
        rC = keypoints(:,3); % Vector to tip of rod 2
        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        % Plot Pole 1
        set(h_pole1,'XData', [0 rB(1)]);
        set(h_pole1,'YData', [0 rB(2)]);

        % Plot Pole 2
        set(h_pole2,'XData' , [rB(1) rC(1)] );
        set(h_pole2,'YData' , [rB(2) rC(2)] );

        pause(.01)
    end
    figure
    plot(tspan,z_out(1,:),'b')
    hold on
    plot(tspan,z_out(2,:),'r')
    xlabel("time [s]")
    ylabel("angle [rad]")
    legend("theta_1","theta_2")
end


function dz = dynamics(z,p)
    % Get mass matrix
    A = A_doublepend(z,p);
    
    % Get forces
    u = [0 0]';
    b = b_doublepend(z,u,p);
    
    % Solve for qdd
    qdd = A\b;
    dz = 0*z;
    
    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end
