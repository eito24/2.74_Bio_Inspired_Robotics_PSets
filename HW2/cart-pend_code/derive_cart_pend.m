name = 'cartpend';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t x dx ddx th dth ddth c l m1 m2 I2 k kappa l0 th0 g F tau real

% Group them
q   = [x; th];      % generalized coordinates
dq  = [dx; dth];    % first time derivatives
ddq = [ddx; ddth];  % second time derivatives
u   = [F; tau];     % controls
p   = [c; l; m1; m2; I2; k; kappa; l0; th0; g];        % parameters

% Generate Vectors and Derivatives
ihat = [0; -1; 0];
jhat = [1; 0; 0];
khat = cross(ihat,jhat);
er2hat =  cos(th)*ihat + sin(th)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rA = x*jhat;
rB = rA+c*er2hat;
rC = rA+l*er2hat;

drA = ddt(rA);
drB = ddt(rB);
drC = ddt(rC);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

T1 = (1/2)*m1*dot(drA, drA);
T2 = (1/2)*m2*dot(drB, drB) + (1/2)* I2 * dth^2;

Vg1 = 0;
Vg2 = m2*g*dot(rB, -ihat);
Ve1 = 1/2*k*(sqrt(simplify(dot(rC,rC)))-l0)^2;
Ve2 = 1/2*kappa*(th - th0)^2;

T = simplify(T1 + T2);
V = Vg1 + Vg2 + Ve1 + Ve2;
Q_tau = M2Q(tau*khat,dth*khat);
Q_F = F2Q(F*jhat,rC);
Q = Q_tau + Q_F;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rA(1:2) rB(1:2) rC(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;


% Rearrange Equations of Motion
A = jacobian(g,ddq);
b = A*ddq - g;

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
