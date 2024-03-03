function zdot = deriv_vel(t,z, coef)
  
    Oxy = partialxy(z);
    
    Ox = Oxy(1);
    Oy = Oxy(2);
    Oz = Oxy(3);
    
    zdot = zeros(6,1);
  
    t_co = coef(5);
  
    beta = coef(1) + coef(2) * (t - t_co) + coef(3) * (t - t_co) ^ 2 + coef(4) * (t - t_co) ^ 3;
    
    DU = 384400;
    TU = 375190;
    mass = 3000; % kg
    thrust = 2; % N
    n0 = (thrust / mass) / 1000;  % km/s^2
    n0 = n0 * TU ^ 2 / DU; 
    cex = 44.645657918834544; % DU/TU
   
    Tm = (cex * n0)/(cex - n0*(t-t_co));
    
    zdot(1) = z(4);
    zdot(2) = z(5);
    zdot(3) = 0;

    temp1 = Tm * cos(beta);
    temp2 = Tm * sin(beta);
    zdot(4) = Ox + 2 * z(5) + temp1;
    zdot(5) = Oy - 2 * z(4) + temp2;
    zdot(6) = 0;
    
end
