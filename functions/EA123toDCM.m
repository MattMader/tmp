function C = EA123toDCM(theta1,theta2,theta3)
% EA123toDCM computes the DCM from the (1-2-3) sequence of Euler angles.
%
% Inputs: 
%   theta1: first Euler angle (1x1) [rad]
%   theta2: second Euler angle (1x1) [rad]
%   theta3: third Euler angle (1x1) [rad]
%
% Outputs:
%   C: computed direction cosine matrix (3x3) [-]
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 21 Feb 2022
%
% Notes:
%   Passed tests against MATLAB's angle2dcm
%

arguments

    theta1 (1,1) {mustBeReal}
    theta2 (1,1) {mustBeReal}
    theta3 (1,1) {mustBeReal}
    
end % arguments

% first rotation (around the 3rd axis with angle theta3)
R3 = [cos(theta3),sin(theta3),0;
    -sin(theta3),cos(theta3),0;
    0,0,1];

% second rotation (around the 2nd axis with angle theta2)
R2 = [cos(theta2),0,-sin(theta2);
    0,1,0;
    sin(theta2),0,cos(theta2)];

% third rotation (around the 1st axis with angle theta1)
R1 = [1,0,0;
    0,cos(theta1),sin(theta1);
    0,-sin(theta1),cos(theta1)];

% final DCM
C = R3*R2*R1;

end % function EA123toDCM