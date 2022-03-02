function theta_dot = KDE_EA313(theta,omega)

A = [sin(theta(3)), cos(theta(3)), 0;
    cos(theta(3))*sin(theta(2)), -sin(theta(3))*sin(theta(2)), 0;
    -sin(theta(3))*cos(theta(2)), -cos(theta(3))*cos(theta(2)), sin(theta(2))];

theta_dot = 1/sin(theta(2))*A*omega;

end