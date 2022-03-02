function C = EPtoDCM(epsilon)
% EPtoDCM computes the DCM from Euler Parameters (quaternion)
%
% Inputs:
%   epsilon: Euler Parameters (4x1) [-]
%
% Outputs:
%   C: Direction cosine matrix (3x3) [-]
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 21 Feb 2022
%
% Notes:
%   

arguments

    epsilon (4,1) {mustBeReal}

end % arguments

% compute DCM elements
C(1,1) = 1-2*epsilon(2)^2-2*epsilon(3)^2;
C(1,2) = 2*(epsilon(1)*epsilon(2)+epsilon(3)*epsilon(4));
C(1,3) = 2*(epsilon(1)*epsilon(3)-epsilon(2)*epsilon(4));
C(2,1) = 2*(epsilon(1)*epsilon(2)-epsilon(3)*epsilon(4));
C(2,2) = 1-2*epsilon(1)^2-2*epsilon(3)^2;
C(2,3) = 2*(epsilon(2)*epsilon(3)+epsilon(1)*epsilon(4));
C(3,1) = 2*(epsilon(1)*epsilon(3)+epsilon(2)*epsilon(4));
C(3,2) = 2*(epsilon(2)*epsilon(3)-epsilon(1)*epsilon(4));
C(3,3) = 1-2*epsilon(1)^2-2*epsilon(2)^2;

end % function EPtoDCM