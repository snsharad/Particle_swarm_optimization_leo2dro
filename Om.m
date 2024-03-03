function omega = Om(x,y,z)
    nu = 0.0121551;
    omega = (x^2 + y^2)/2 + (1 - nu)/sqrt((x + nu)^2 + y^2 + z^2) + nu/sqrt((x + nu -1)^2 + y^2 + z^2);
end