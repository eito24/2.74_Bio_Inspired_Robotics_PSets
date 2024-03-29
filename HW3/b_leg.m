function b = b_leg(in1,in2,in3)
%B_LEG
%    B = B_LEG(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    02-Oct-2022 21:44:53

dth1 = in1(3,:);
dth2 = in1(4,:);
g = in3(19,:);
l_AC = in3(17,:);
l_A_m3 = in3(13,:);
l_B_m2 = in3(12,:);
l_C_m4 = in3(14,:);
l_OA = in3(15,:);
l_OB = in3(16,:);
l_O_m1 = in3(11,:);
m1 = in3(1,:);
m2 = in3(2,:);
m3 = in3(3,:);
m4 = in3(4,:);
tau1 = in2(1,:);
tau2 = in2(2,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = sin(th1);
t3 = sin(th2);
t4 = th1+th2;
t5 = dth1.^2;
t6 = l_OA.*t2;
t7 = sin(t4);
mt1 = [tau1+tau2+dth2.*(dth1.*l_AC.*l_C_m4.*m4.*t3.*2.0+dth2.*l_AC.*l_C_m4.*m4.*t3+dth1.*l_AC.*l_OA.*m4.*t3.*2.0+dth2.*l_AC.*l_OA.*m4.*t3+dth1.*l_A_m3.*l_OA.*m3.*t3.*2.0+dth2.*l_A_m3.*l_OA.*m3.*t3+dth1.*l_B_m2.*l_OB.*m2.*t3.*2.0+dth2.*l_B_m2.*l_OB.*m2.*t3)-g.*m2.*(l_B_m2.*t7+l_OB.*t2)-g.*m3.*(t6+l_A_m3.*t7)-g.*m4.*(t6+l_AC.*t7+l_C_m4.*t2)-g.*l_O_m1.*m1.*t2];
mt2 = [tau2+dth2.*(dth1.*l_AC.*l_C_m4.*m4.*t3+dth1.*l_AC.*l_OA.*m4.*t3+dth1.*l_A_m3.*l_OA.*m3.*t3+dth1.*l_B_m2.*l_OB.*m2.*t3)-g.*l_AC.*m4.*t7-g.*l_A_m3.*m3.*t7-g.*l_B_m2.*m2.*t7-l_AC.*l_C_m4.*m4.*t3.*t5-l_AC.*l_OA.*m4.*t3.*t5-l_A_m3.*l_OA.*m3.*t3.*t5-l_B_m2.*l_OB.*m2.*t3.*t5-dth1.*dth2.*l_AC.*l_C_m4.*m4.*t3-dth1.*dth2.*l_AC.*l_OA.*m4.*t3-dth1.*dth2.*l_A_m3.*l_OA.*m3.*t3-dth1.*dth2.*l_B_m2.*l_OB.*m2.*t3];
b = [mt1;mt2];
