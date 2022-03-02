function [theta1,theta2,theta3] = DCMtoEA323(C)
% DCMtoEA323 computes the (3-2-3) sequence of Euler angles from the DCM.
%
% Inputs:
%   C: provided direction cosine matrix (3x3) [-]
%
% Outputs:
%   theta1: first Euler angle (1x1) [rad]
%   theta2: second Euler angle (1x1) [rad]
%   theta3: third Euler angle (1x1) [rad]
%

arguments
    
    C (3,3) {mustBeReal}

end % arguements

% first Euler angle
theta1 = atan2(C(3,2),C(3,1)); % [rad]

% second Euler angle
theta2 = acos(C(3,3)); % [rad]

% third Euler angle
theta3 = atan2(C(2,3),-C(1,3)); % [rad]

end % function DCMtoEA323