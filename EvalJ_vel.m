function J = EvalJ_vel(P, J, N_particles)
      
    parfor it = 1:N_particles 

        c0 = P(it,1);
        c1 = P(it,2);
        c2 = P(it,3);
        c3 = P(it,4);
        
        t_spi = 38.95;
        t_co = P(it,5);
        t_f = P(it,6);
        
        c = [c0 c1 c2 c3 t_co];
        
        if t_f <= t_co
            J(it) = 1000;
            continue;
        end
        
        options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

        %%%%%%%% Escape spiral

        % Corresponds to:
        %         const1 = 375190 / 384400;
        %         rp = 6378 + 350;
        %         e = 0;
        %         a = rp / (1 - e);
        %         i = 0;
        %         argPeriap = 0;
        %         raan = 0;
        %         trueAnom = pi/12;
                
        %       mee_initial = oe2mee(a, e, i, argPeriap, raan, trueAnom); 
        %         tspan = [0 t_spi];
        %         [~, mee] = ode45(@meeEoms, tspan, mee_initial, options);  
        %         eci_init_co = mee(end,:);
       
        %%%%%%%%%% Coasting arc
        
        eci_init_co = [-0.418755220934786  0.670814669276861  0  -0.118988146399635  1.230920167295980  0]';

        tspan = [t_spi t_co];
        [~, eci_co] = ode45(@deriv, tspan, eci_init_co, options);
        
        eci_init_trf = eci_co(end,:);
                        
        tspan = [t_co t_f];
        [~, eci_trf] = ode45(@(t, mee_trf)deriv_vel(t, mee_trf, c), tspan, eci_init_trf, options);
                            
        xT = eci_trf(end,1);
        yT = eci_trf(end,2);
        vxT = eci_trf(end,4);
        vyT = eci_trf(end,5);  

        
        % Target DRO data point
        xF = 1.678503;
        yF = 0.129033;
        vxF = 0.1232861;
        vyF = -1.255867;
         
        d01 = vxF - vxT;
        d02 = vyF - vyT;
        d1 = abs(xF - xT);
        d2 = abs(yF - yT);
        d4 = abs(vxF - vxT);
        d5 = abs(vyF - vyT);

        tol = 1.0e-12;

        if abs(d1) > tol
            a1 = 1000;
        else 
            a1 = 1;
        end
        
        if abs(d2) > tol
            a2 = 1000;
        else 
            a2 = 1;
        end
        
        if abs(d4) > tol
            a4 = 100;
        else 
            a4 = 1;
        end
        
        if abs(d5) > tol
            a5 = 100;
        else 
            a5 = 1;
        end

       meanError = a1*(d1^2) + a2*(d2^2) + a4*(d4^2) + a5*(d5^2);
       meanError = meanError / 4;
       J(it) = sqrt(meanError);

        if  it == N_particles
            fprintf('\nEvalJ: % d %e\n %e  %e  %e  %e\n',it, J(it), d1, d2, d01, d02);
        end
    end     
end



