function A = A_cartpend(in1,in2)
%A_cartpend
%    A = A_cartpend(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    27-Sep-2022 22:36:46

I2 = in2(5,:);
c = in2(1,:);
m1 = in2(3,:);
m2 = in2(4,:);
th = in1(2,:);
t2 = cos(th);
t3 = c.*m2.*t2;
A = reshape([m1+m2,t3,t3,I2+c.^2.*m2],[2,2]);