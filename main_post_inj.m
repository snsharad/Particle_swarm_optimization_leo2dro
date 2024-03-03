

X = csvread('GBest_xfer.csv');

rp = 6378 + 500;
e = 0;
a = rp / (1 - e);
i = 0;
argPeriap = 0;
raan = 0;
trueAnom = 0;
%c = [0.001 1.0e-08 -5e-10 -5e-08 ];


tu = 375190;
tspi = 38.95;
tco = X(5);
tf = X(6);
tj = tf + 6;


% c = [-0.000097627936475  -0.000803127174647  -0.000144217005735  -0.00057846239691  tco];

c = [X(1)  X(2)  X(3)  X(4) tco];

mee_initial = oe2mee(a, e, i, argPeriap, raan, trueAnom);

tspan = [0 tspi];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t, mee] = ode45(@meeEoms, tspan, mee_initial, options);


ht = length(mee);

eci = zeros(ht,6);

angle = zeros(ht,1);

% Convert from modified equinotial elements to
%	earth centric x, y, z, xdot, ydot, zdot
% spiral
for j = 1:ht

    eci(j,:) = mee2eci(mee(j,:));
   
    f = mee(j,2);
    g = mee(j,3);
    L = mee(j,6);
    e = sqrt(f ^ 2 + g ^ 2);
    w = 1 + f * cos(L) + g * sin(L);     
    temp1 = (e ^ 2 + 2 * e - w ^ 2 + 1);
    temp2 = w * (f * sin(L) - g * cos(L));
    angle(j) = atan(temp1 / temp2);
    
end



% Angular velocity
omega = [0 0 1.013013]; % in rad/TU

z = zeros(ht,3);
v = zeros(ht,3);
z1 = zeros(ht,3);
v1 = zeros(ht,3);
 
 % Non-dimensional
 % Change to rotating coordinates
 % spiral
 for q = 1:ht
   
   %t(q) = t(q) / 375190;
   z(q, 1) = eci(q,1) / 384400; % DU
   z(q, 2) = eci(q,2) / 384400; % DU
   z(q, 3) = 0;
   
   v(q, 1) = eci(q,4) * 375190 / 384400;
   v(q, 2) = eci(q,5) * 375190 / 384400;
   v(q, 3) = 0;
   
   w = omega(3) * t(q);  
   
   C_IB = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];  
   
   int = C_IB * z(q,:)';
%    int_v = C_IB * v(q,:)';
   cr_prod = cross(omega, z(q,:)); 
   int_v = v(q,:) - cr_prod; % convert vel magnitude to rotating frame
   int_v = C_IB * int_v';    % multiply by DCM to get the correct direction   
   
   z1(q,:) = [int(1) int(2) int(3)];
   v1(q,:) = [int_v(1) int_v(2) int_v(3)];
 end
 
x1 = [(z1(:,1)-0.0121551) z1(:,2) v1(:,1) v1(:,2)];

coast_init = [x1(end,1); x1(end,2); 0; x1(end,3); x1(end,4); 0];
tspan = [tspi tco];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t_co, eci_co] = ode45(@deriv, tspan, coast_init, options);
 
trf_init = eci_co(end,:)';

tspan = [tco tf];
[t_trf, eci_trf] = ode45(@(t_trf, eci_trf)deriv_vel(t_trf, eci_trf, c), tspan, trf_init, options);

post_inj_init = eci_trf(end,:)';

tspan = [tf tj];
[t_inj, inj] = ode45(@deriv, tspan, post_inj_init, options);

ht1 = length(t_trf);
angle = zeros(ht1,1);


for iter = 1:ht1
    angle(iter) = c(1) + c(2) * (t_trf(iter) - tco) + c(3) * (t_trf(iter) - tco) ^ 2 + c(4) * (t_trf(iter) - tco) ^ 3;
    angle(iter) = angle(iter) * 180 / pi;
end



%%
%%
moon = 1 - 0.0121551 ;
L1 = 0.836893; 
L2 = 1.14878370300839;
N = csvread('DRO_2.csv');
%%
figure(1)
plot(x1(:,1), x1(:,2));
hold on;
plot(x1(1,1), x1(1,2), '*');
plot(x1(end,1), x1(end,2), '+');
plot(eci_co(:,1), eci_co(:,2));
plot(eci_trf(1,1), eci_trf(1,2), '*');
plot(eci_trf(:,1), eci_trf(:,2), 'green');
plot(eci_trf(end,1), eci_trf(end,2), '+');
plot(inj(:,1), inj(:,2), 'red');
plot(inj(end,1), inj(end,2), '+');

plot(N(:,1), N(:,2));
plot(moon, 0, 'o');
plot(-0.0121552, 0, 'o');
plot(L1, 0, '*');
plot(L2, 0, '*');
xlabel('X (DU)');
ylabel('Y (DU)');
grid on;
axis equal;
hold off;

% figure(2)
% plot(eci(:,1), eci(:,2));
% 
figure(3)
plot(t_trf(:), angle(:));
grid on;
xlabel('Thrusting time (TU)');
ylabel('Angles (deg)');

% figure(2)
% plot(t(:), angle(:));
