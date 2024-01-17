function foot_position = position_foot(in1,in2)
%POSITION_FOOT
%    FOOT_POSITION = POSITION_FOOT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    02-Oct-2022 21:44:54

l_AC = in2(17,:);
l_DE = in2(18,:);
l_OB = in2(16,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
foot_position = [l_DE.*t3+l_OB.*t3+l_AC.*sin(t4);-l_DE.*t2-l_OB.*t2-l_AC.*cos(t4);0.0];
