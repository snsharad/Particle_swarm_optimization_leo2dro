
clc;
clear;

format long;
warning('off', 'all');
toli = 1.0e-12;    
DU = 384400; 
         
% Define constants
m2 = 7.346*1e22 ;
m1 = 5.9724 *1e24;
mu = m2/(m1+m2);

L1 = 0.836915727639612; 
L2 = 1.155681694999380;
moon = 1 - mu ;

%% PSO Basic

N_particles = 100;
 
%% Number of variables
N_elements = 6;
N_iterations = 1000;

% For reference:
% Variable 1 : c0 
% Variable 2 : c1 
% Variable 3 : c2
% Variable 4 : c3
% Variable 5 : t_co
% Variable 6 : t_f

% Set lower and upper bounds on unknowns (particle elements)

BLp = [-3, -3, -3, -3, 39.5, 41];

% Lower and upper bound for T, eqn 31
BUp = [3, 1, 1, 0, 40.5, 42.5];

%%
P = zeros(N_particles,N_elements);
    
%     Create random initial population, with each particle between BLp and BUp
   for i=1:N_elements
      P(:,i) = BLp(i) + rand(N_particles,1)*(BUp(i)-BLp(i));
   end
    
%     Write initial swarm to spreadsheet
    filename = 'Initial_swarm_xfer.csv';
%     
    dlmwrite(filename, P, 'precision', '%e');
    
    PBest = zeros(N_particles,N_elements);
    J = zeros(N_particles,1);  % values of J for the entire population
    JBest = zeros(N_particles,1);  %  values of J corresponding to PBest values
    V = zeros(N_particles, N_elements);  
    GGstar = zeros(1,N_iterations);
    w = zeros(N_iterations, 2);
    
    %% determine velocity bounds
    BUv = BUp - BLp;
    BLv = -BUv;
    
    for i = 1:N_particles
        JBest(i) = inf;
    end
    
    GG = inf;   % value of J for the global best particle GBest
  
    for j = 1:1:N_iterations
       
        J = EvalJ_vel(P,J,N_particles);
        
        %% Update PBest (for each particle) and GBest (for the population)
        for i = 1:N_particles
            if J(i) < JBest(i)
                PBest(i,:) = P(i,:);
                JBest(i) = J(i);
            end
            if J(i) < GG
                GG = J(i);
                GBest = P(i,:);
            end
        end
        
        %% Update velocity for each particle
        c_I = (1 + rand)/2;
        c_C = 1.49445*rand;
        c_S = 1.49445*rand;
        for i =1:N_particles
            V(i,:) = c_I*V(i,:) + c_C*(PBest(i,:) - P(i,:)) + c_S*(GBest - P(i,:));
            for k = 1:N_elements
                if V(i,k) < BLv(k)
                    V(i,k) = BLv(k);
                elseif V(i,k) > BUv(k)
                    V(i,k) = BUv(k);
                end
            end
        end
        
        %% Update position for each particle
        for i =1:N_particles
            P(i,:) = P(i,:) + V(i,:);
            for k = 1:N_elements
                if P(i,k) < BLp(k)
                    P(i,k) = BLp(k);
                    V(i,k) = 0;
                elseif P(i,k) > BUp(k)
                    P(i,k) = BUp(k);
                    V(i,k) = 0;
                end
            end
        end
      
        GGstar(j) = GG; % global best J after each iteration j
%        if mod(j,20)==0
             fprintf('\n\n%d  %e\n', j, GG);   
             fprintf('%e  ', GBest); 
%        end 
       
        w(j,:) = [j, GG];
    end
    
    fprintf('\n');
    fprintf('%e  ', GBest);
             
    semilogy(w(:,1), w(:,2)); % Plot J vs iteration
  
    filename = 'GBest_xfer.csv';

    dlmwrite(filename, GBest, 'precision', '%e');

