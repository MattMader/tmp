function epsilon_dot = KDE_EP(epsilon,omega)

A = [0,omega(3),-omega(2),omega(1)
    -omega(3),0,omega(1),omega(2)
    omega(2),-omega(1),0,omega(3)
    -omega(1),-omega(2),-omega(3),0];

epsilon_dot = 0.5*A*epsilon;

end % function