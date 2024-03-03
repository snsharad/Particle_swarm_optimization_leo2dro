function eci = mee2eci(m)
    
    mu = 398600.6618;
    
    p = m(1);
    f = m(2);
    g = m(3);
    h = m(4);
    k = m(5);
    L = m(6);
    
    alpha_sq = h^2 - k^2;
    s_sq = 1 + h^2 + k^2;
    w = 1 + f * cos(L) + g * sin(L);
    r = p / w;
    
%     eci = zeros(6,1);
    
    temp1 = r / s_sq;
    temp2 = (1 / s_sq) * sqrt(mu / p);
    
    eci(1) = temp1 * (cos(L) + alpha_sq * cos(L) + 2 * h * k * sin(L));
    eci(2) = temp1 * (sin(L) - alpha_sq * sin(L) + 2 * h * k * cos(L));
    eci(3) = 2 * temp1 * (h * sin(L) - k * cos(L));
    eci(4) = - temp2 * (sin(L) + alpha_sq * sin(L) - 2 * h * k * cos(L) + g - 2 * f * h * k + alpha_sq * g);
    eci(5) = - temp2 * (- cos(L) + alpha_sq * cos(L) + 2 * h * k * sin(L) - f + 2 * g * h * k + alpha_sq * f);
    eci(6) = 2 * temp2 * (h * cos(L) + k * sin(L) + f * h + g * k);
    
end