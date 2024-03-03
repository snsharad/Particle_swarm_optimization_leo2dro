function mee = oe2mee(a, e, i, argPeriap, raan, trueAnom) 

    mee = zeros(6,1);
    
    mee(1) = a * (1 - e^2);                 %semi latus rectum
    mee(2) = e * cos(argPeriap + raan);     % f
    mee(3) = e * sin(argPeriap + raan);     % g
    mee(4) = tan(i/2) * cos(raan);          % h
    mee(5) = tan(i/2) * sin(raan);          % k
    mee(6) = raan + argPeriap + trueAnom;   % L
    
end