function zdot = deriv(t,z)

    Oxy = partialxy(z);
    
    Ox = Oxy(1);
    Oy = Oxy(2);
    Oz = Oxy(3);
    
    zdot = zeros(6,1);
    
    zdot(1) = z(4);
    zdot(2) = z(5);
    zdot(3) = z(6);
    zdot(4) = Ox + 2 * z(5);
    zdot(5) = Oy - 2 * z(4);
    zdot(6) = Oz;
     
end