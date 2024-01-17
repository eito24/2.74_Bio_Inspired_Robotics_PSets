name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t th1 dth1 ddth1 th2 dth2 ddth2 m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g tau1 tau2 real

% Group them
q   = [th1; th2];      % generalized coordinates
dq  = [dth1; dth2];    % first time derivatives
ddq = [ddth1; ddth2];  % second time derivatives
u   = [tau1; tau2];     % controls
p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';       % parameters

% Generate Vectors and Derivatives
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);
erBhat = sin(th1)*ihat - cos(th1)*jhat;    %direction that body 1 and body 3 extends
erChat = sin(th1+th2)*ihat - cos(th1+th2)*jhat; %direction that body 2 extends

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rA = l_OA*erBhat;
rB = l_OB*erBhat;
rC = rA + l_AC*erChat;
rD = rB + l_AC*erChat;
rE = rD + l_DE*erBhat;
r1c = l_O_m1*erBhat; % center of mass for body 1
r2c = rB + l_B_m2*erChat;
r3c = rA + l_A_m3*erChat;
r4c = rC + l_C_m4*erBhat;


drA = ddt(rA);
drB = ddt(rB);
drC = ddt(rC);
drD = ddt(rD);
drE = ddt(rE);
dr1c = ddt(r1c);
dr2c = ddt(r2c);
dr3c = ddt(r3c);
dr4c = ddt(r4c);


% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces

T1 = (1/2) * m1 * dot(dr1c, dr1c) + (1/2)* (I1) * dth1^2; % KE of body 1
T2 = (1/2) * m2 * dot(dr2c, dr2c) + (1/2)* (I2) * (dth1+dth2)^2; % KE of body 2
T3 = (1/2) * m3 * dot(dr3c, dr3c) + (1/2)* I3 * (dth1+dth2)^2; % KE of body 3
T4 = (1/2) * m4 * dot(dr4c, dr4c) + (1/2)* I4 * dth1^2; % KE of body 4
T5 = (1/2) * Ir * (dth1*N)^2; % KE of Rotor 1
T6 = (1/2) * Ir * (dth1+N*dth2)^2; % KE of Rotor 2

Vg1 = m1*g*dot(r1c, jhat); % GPE body 1
Vg2 = m2*g*dot(r2c, jhat);
Vg3 = m3*g*dot(r3c, jhat);
Vg4 = m4*g*dot(r4c, jhat);


T = simplify(T1 + T2 + T3 + T4 + T5 + T6);
V = Vg1 + Vg2 + Vg3 + Vg4;

Q_tau1 = M2Q(tau1*khat,dth1*khat);
Q_tau2 = M2Q(tau2*khat,(dth1+dth2)*khat);

%Q_F = F2Q(F*jhat,rC);
Q = Q_tau1+Q_tau2;
%Q = [tau1;tau2;];

% Assemble the array of cartesian coordinates of the key points
keypoints = [rA(1:2) rB(1:2) rC(1:2) rD(1:2) rE(1:2) r1c(1:2) r2c(1:2) r3c(1:2) r4c(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+V;
L = T-V;
G = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

% Rearrange Equations of Motion
A = jacobian(G,ddq);
b = A*ddq - G;
foot_position = rE;
jacobian_rE = jacobian(rE,q);

% Write Energy Function and Equations of Motion
z  = [q ; dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(T, 'file', ['T_' name], 'vars', {z, p});
matlabFunction(V, 'file', ['V_' name], 'vars', {z, p});
matlabFunction(jacobian_rE, 'file', ['jacobian_rE'], 'vars', {z, p});
matlabFunction(foot_position,'file','position_foot','vars',{z p});
matlabFunction(drE,'file','drE','vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});
