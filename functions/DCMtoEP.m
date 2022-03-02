function epsilon = DCMtoEP(C)
% DCMtoEP calculates the Euler Parameters (quaternion) from a DCM using
% [Sheppard 1978] method for numerical stability. Positive epsilon4
%
% Inputs:
%   C: Direction cosine matrix (3x3) [-]
%
% Outputs:
%   epsilon: Euler Parameters (4x1) [-]
%
% Information:
%   Author: Matthew Mader
%   Contact: maderm@purdue.edu
%   Date: 21 Feb 2022
%
% Notes:
%   Passed tests against MATLAB's dcm2quat function.
%

arguments

    C (3,3) {mustBeReal}
    
end % arguments

% initialize output
epsilon = zeros(4,1);

% find largest value
[v,i] = max( ...
    [0.25*(1+2*C(1,1)-trace(C)); ...
    0.25*(1+2*C(2,2)-trace(C));
    0.25*(1+2*C(3,3)-trace(C));
    0.25*(1+trace(C))]);

% largest value
switch i
    % epsilon1^2
    case 1
        epsilon(1) = sqrt(v); % [-]
        epsilon(2) = (C(1,2)+C(2,1))/(4*epsilon(1)); % [-]
        epsilon(3) = (C(3,1)+C(1,3))/(4*epsilon(1)); % [-]
        epsilon(4) = (C(2,3)-C(3,2))/(4*epsilon(1)); % [-]
        
    % epsilon2^2
    case 2
        epsilon(2) = sqrt(v); % [-]
        epsilon(1) = (C(1,2)+C(2,1))/(4*epsilon(2)); % [-]
        epsilon(3) = (C(2,3)+C(3,2))/(4*epsilon(2)); % [-]
        epsilon(4) = (C(3,1)-C(1,3))/(4*epsilon(2)); % [-]
        
    % epsilon3^2
    case 3
        epsilon(3) = sqrt(v); % [-]
        epsilon(1) = (C(3,1)+C(1,3))/(4*epsilon(3)); % [-]
        epsilon(2) = (C(2,3)+C(3,2))/(4*epsilon(3)); % [-]
        epsilon(4) = (C(1,2)-C(2,1))/(4*epsilon(3)); % [-]
        
    % epsilon4^2
    otherwise
        epsilon(4) = sqrt(v); % [-]
        epsilon(1) = (C(2,3)-C(3,2))/(4*epsilon(4)); % [-]
        epsilon(2) = (C(3,1)-C(1,3))/(4*epsilon(4)); % [-]
        epsilon(3) = (C(1,2)-C(2,1))/(4*epsilon(4)); % [-]

end % switch


% if constraint variable, epsilon4, is negative
if epsilon(4) < 0

    epsilon = -epsilon;

end % if


end % function DCMtoEP