function mdot = meeEoms(t, m)
    
    mu = 398600.6618; % gravitational param of Earth
    p = m(1);
    f = m(2);
    g = m(3);
    h = m(4);
    k = m(5);
    L = m(6);
    
    tu = 375190;
    xy = mee2eci(m);   
    sm = [384400 0 0]; 
    
    rad = xy(1:3);
    vel = xy(4:6);
    
    d = rad - sm;
    dMag = norm(d);
    
    smMag = norm(sm);
      
    term1 = d / dMag ^3 + sm / smMag ^3;
    moonPert = - 4902.8695 * term1; % mu of moon
    
    ihat_r = rad / norm(rad);
    
    term2 = cross(rad, vel);
    ihat_n = term2 / norm(term2);
    
    ihat_t = cross(ihat_n, ihat_r);
    
    Q = [ihat_r ihat_t ihat_n];
 
    pert = Q' * moonPert;
 
    
    mdot = zeros(6,1);
   
    e = sqrt(f^2 + g^2);
    s_squared = 1 + h^2 + k^2;
    w = 1 + f * cos(L) + g * sin(L);
    r = p / w;
    
    tempA = e^2 + 2 * e - w^2 + 1;
    tempB = w * (f * sin(L) - g * cos(L));
    
    if (tempA ~= 0) || (tempB ~= 0)
        beta = atan(tempA / tempB);
    else
        beta = 0;
    end
    
    mass = 3000; % kg
    thrust = 2; % N
    n0 = (thrust / mass) / 1000; % km/s^2
    c = 45.7416;
    Tm = (c * n0)/(c - n0*t*tu);
      
    acc_r = Tm * sin(beta) + pert(1);
    acc_t = Tm * cos(beta) + pert(2);
    acc_n = 0.0;
        
    temp1 = sqrt(p / mu);
    temp2 = (w + 1) * cos(L) + f; 
    temp3 = (w + 1) * sin(L) + g;
    temp4 = h * sin(L) - k * cos(L);
    temp5 = s_squared / (2 * w);
    
    tmp = temp1 * 2 * r * acc_t;
    mdot(1) = tu * tmp;
    tmp = temp1 * (acc_r * sin(L) + temp2 * acc_t / w - temp4 * g * acc_n / w);
    mdot(2) = tu * tmp;
    tmp = temp1 * (- acc_r * cos(L) + temp3 * acc_t / w + temp4 * f * acc_n / w);
    mdot(3) = tu * tmp;
    
    mdot(4) = tu * temp1 * temp5 * acc_n * cos(L);
    mdot(5) = tu * temp1 * temp5 * acc_n * sin(L);
    tmp = sqrt(mu * p) * (1 / r)^2 + (1 / w) * temp1 * temp4 * acc_n;
    mdot(6) = tu * tmp;
    
end