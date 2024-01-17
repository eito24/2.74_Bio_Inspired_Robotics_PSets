function simulate_cartpend()

    %% Definte fixed paramters
    c= .5;
    l = 1;
    m1 = 1;
    m2 = 1;
    I2 = .1;
    k = 3;
    kappa = .5;
    l0 = 0;
    th0 = 0;
    g = 10;

    p   = [c; l; m1; m2; I2; k; kappa; l0; th0; g];        % parameters

    %% Perform Dynamic simulation    
    dt = 0.001;
    tf = 15;
    num_steps = floor(tf/dt);
    tspan = linspace(0, tf, num_steps); 
    z0 = [1; 0; 0; 0];
    z_out = zeros(4,num_steps);
    z_out(:,1) = z0;
    for i=1:num_steps-1
        dz = dynamics(z_out(:,i), p);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
%         z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;
%         z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt; % + 0.5*dz(3:4)*dt*dt;
    end
    final_state = z_out(:,end);
    

    %% Compute Energy
    E = energy_cartpend(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');

    %% Animate Solution
    figure(2); clf;
        % Prepare plot handles
    hold on
    h_cart = plot([0],[0],'LineWidth',5);
    h_pole = plot([0],[0],'LineWidth',3);
    h_spring = plot([0],[0],'k');
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
        keypoints = keypoints_cartpend(z,p);

        rA = keypoints(:,1); % Vector to base of cart
        rC = keypoints(:,3); % Vector to tip of pendulum

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        % Plot Cart
        set(h_cart,'XData', [ rA(1) - .2 ; rA(1) + .2]);
        set(h_cart,'YData', [0 0]);

        % Plot Pole
        set(h_pole,'XData' , [rA(1) rC(1)] );
        set(h_pole,'YData' , [rA(2) rC(2)] );

        % Plot spring
        set(h_spring, 'XData' , [0 rC(1)] );
        set(h_spring, 'YData' , [0 rC(2)] );

        pause(.01)
    end
end

function dz = dynamics(z,p)
    % Get mass matrix
    A = A_cartpend(z,p);
    
    % Get forces
    u = [0 0]';
    %u = [ 0 -5*z(4)]';
    b = b_cartpend(z,u, p);
    
    % Solve for qdd
    qdd = A\b;
    dz = 0*z;
    
    % Form dz
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end
