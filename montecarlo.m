"matlab.mlintpath"

% One Dimensional Brownian Motion
N = 1000;
displacement = randn(1,N);
plot(displacement);

% Distribution of Displacements
hist(displacement, 25);

% Convert displacements to position
x = cumsum(displacement);
plot(x);
ylabel('position');
xlabel('time step');
title('Position of 1D Particle versus Time');

% Two Dimensional Particle Simulation
particle = struct();
particle.x = cumsum( randn(N, 1) );
particle.y = cumsum( randn(N, 1) );
plot(particle.x, particle.y);
ylabel('Y Position');
xlabel('X Position');
title('position versus time in 2D');

% Compute the Displacement Squared
dsquared = particle.x .^ 2 + particle.y .^ 2;
plot(dsquared);

% Theoretical Value of D
d    = 1.0e-6;              % diameter in meters
eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
kB   = 1.38e-23;            % Boltzmann constant
T    = 293;                 % Temperature in degrees Kelvin

D    = kB * T / (3 * pi * eta * d);

%ans  = 4.2902e-013;

% A more realistic particle 
dimensions = 2;         % two dimensional simulation
tau = .1;               % time interval in seconds
time = tau * (1:N);       % create a time vector for plotting

k = sqrt(D * dimensions * tau);
dx = k * randn(N,1);
dy = k * randn(N,1);

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
 squaredDisplacement = ( x .^ 2) + ( y .^ 2);

plot(x,y);
title('Particle Track of a Single Simulated Particle');

% Displacement Squared Plot
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line

plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time for 1 Particle in 2 Dimensions');

% Estimating D from the Simulated Data
simulatedD = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau );

%ans = 4.2192e-013;

% Uncertainty in the Estimate
standardError = std( dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) );
actualError = D - simulatedD;

%standardError = 1.3162e-014

%actualError   = 7.1019e-015

% Systematic Error -- Bulk Flow in the Solvent
dx = dx + 0.2 * k;
dy = dy + 0.05 * k;

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
squaredDisplacement = ( x .^ 2) + ( y .^ 2);

simulatedD    = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau )
standardError = std(  dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) )
actualError = D - simulatedD

plot(x,y);
title('Particle Track of a Single Simulated Particle with Bulk Flow');

%simulatedD    =  4.2926e-013
%standardError =  1.3694e-014
%actualError   =  -2.3859e-016

% Displacement Squared in the Presence of Bulk Flow
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line
plot(time, squaredDisplacement);
hold off;

xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time with Bulk Flow');

% Simulating Multiple Particles

% Generating Multiple Data Sets
particleCount = 10;
N = 50;
tau = .1;
time = 0:tau:(N-1) * tau;
particle = { };             % create an empty cell array to hold the results

for i = 1:particleCount
    particle{i} = struct();
    particle{i}.dx = k * randn(1,N);
    particle{i}.x = cumsum(particle{i}.dx);
    particle{i}.dy = k * randn(1,N);
    particle{i}.y = cumsum(particle{i}.dy);
    particle{i}.drsquared = particle{i}.dx .^2 + particle{i}.dy .^ 2;
    particle{i}.rsquared = particle{i}.x .^ 2 + particle{i}.y .^ 2;
    particle{i}.D = mean( particle{i}.drsquared ) / ( 2 * dimensions * tau );
    particle{i}.standardError = std( particle{i}.drsquared ) / ( 2 * dimensions * tau * sqrt(N) );
end

% Results
clf;
particleCount = 10;
hold on;
for i = 1:particleCount
    plot(particle{i}.x, particle{i}.y, 'color', rand(1,3));
end

xlabel('X position (m)');
ylabel('Y position (m)');
title('Combined Particle Tracks');
hold off;

% Displacement Squared
% compute the ensemble average
rsquaredSum = zeros(1,N);

for i = 1:particleCount
    rsquaredSum = rsquaredSum + particle{i}.rsquared;
end

ensembleAverage = rsquaredSum / particleCount;

% create the plot
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'b', 'LineWidth', 3);      % plot theoretical line

plot(time, ensembleAverage , 'k', 'LineWidth', 3);          % plot ensemble average
legend('Theoretical','Average','location','NorthWest');

for i = 1:particleCount
    plot(time, particle{i}.rsquared, 'color', rand(1,3));   % plot each particle track
end

xlabel('Time (seconds)');
ylabel('Displacement Squared (m^2)');
title('Displacement Squared vs Time');
hold off;

% Estimated Value of D
clear D e dx;

% extract the D value from each simulation and place them all into a single
% matrix called 'D'
for i = 1:particleCount
    D(i) = particle{i}.D;
    dx(i,:) = particle{i}.dx;
    e(i) = particle{i}.standardError;
end

% compute the estimate of D and the uncertainty
averageD = mean(D)
uncertainty = std(D)/sqrt(particleCount)

% plot everything
clf;
hold on;

plot(averageD * ones(1,particleCount), 'b', 'linewidth', 3);                    % plot estimated D
plot((averageD + uncertainty) * ones(1,particleCount), 'g-', 'linewidth', 1);   % plot upper error bar
plot((averageD - uncertainty) * ones(1,particleCount), 'g-', 'linewidth', 1);   % plot lower error bar
errorbar(D,e,'ro');                                                             % plot D values with error bars

xlabel('Simulation Number');
ylabel('Estimated Diffusion Coefficient');
title('Estimated Diffusion Coefficient with Error Bars')
legend('Average Value of D', 'location', 'NorthWest');

hold off;

%averageD    = 4.2886e-013

%uncertainty = 2.3294e-014

% More Advanced Statistics and Plots
% A More Perfect Distribution
hist(randn(1,1e3),25)
xlabel('Value');
ylabel('Frequency');
title('Histogram of Values for 1000 Samples of randn');

hist(randn(1,1e6),100)
xlabel('Value');
ylabel('Frequency');
title('Histogram of Values for 1000000 Samples of randn');

% Sampling Uncertainty
clf;
hold on;
<span class="keyword">for i= 1:100
    <span class="keyword">for j = 1:50
        y = randn(1,i);
        m = mean(y);
        plot(i,m,'x');
    <span class="keyword">end

<span class="keyword">end
plot(1:100, 1./sqrt(1:100), 'k', 'LineWidth', 2);   % plot upper error bar in dark black
plot(1:100, -1./sqrt(1:100), 'k', 'LineWidth', 2);  % plot lower error bar in dark black
hold off;
xlabel('Population Size (N)');
ylabel('Population Mean');
title('Population Mean of Normally Distributed Random Variable versus Population Size');

% Squaring a Random Variable?
dx = randn(1,1e6);
dy = randn(1,1e6);
drSquared = dx .^ 2 + dy .^ 2;

mean(drSquared)
var(dx) + var(dy)

clf;
hist(drSquared,100);
title('Chi Squared Distribution (2 DOF), 1000000 Samples');

%ans = 1.9982
%ans = 1.9982

% 100 Years of BIO Lab Data in 1 Second
clf;
hold on;

for i= 1:100
    for j = 1:50
        dx = randn(1,i);
        dy = randn(1,i);
        m = mean( dx .^ 2 + dy .^ 2 );
        plot(i,m,'x');
    end
end

plot(1:100, 2 + 2./sqrt(1:100), 'k', 'LineWidth', 2);  % plot upper error bar in dark black
plot(1:100, 2 - 2./sqrt(1:100), 'k', 'LineWidth', 2);  % plot lower error bar in dark black

hold off;
xlabel('Population Size (N)');
ylabel('Population Mean');
title('Population Mean of Chi Squared (2 DOF) Distributed Random Variable versus Population Size');

% Auto Correlation 
clf;
c = xcorr(particle{1}.dx, 'coeff');
xaxis = (1-length(c))/2:1:(length(c)-1)/2;

plot(xaxis, c);
xlabel('Lag');
ylabel('Correlation');
title('Particle 1 x-axis Displacement Autocorrelation');

% Cross Correlation
clf;
c = xcorr(particle{1}.dx, particle{2}.dx, 'coeff');
xaxis = (1-length(c))/2:1:(length(c)-1)/2;
plot(xaxis, c);
xlabel('Lag');
ylabel('Correlation');
Title('Particle 1, 2 x-axis Displacement Cross Correlation');

% Advanced Methods
% create an array whose columns contain the dx values for each particle

for i = 1:particleCount
    allDx(:,i) = particle{i}.dx';
end

% compute all possible auto and cross correlations

c = xcorr(allDx, 'coeff');

% plot the results
clf;
hold on;
for i=1:size(c,1)
    plot(xaxis, c(:,i),'color',rand(1,3));
end
hold off;

xlabel('Lag');
ylabel('Correlation Coefficient');
title('All Possible Auto and Cross Correlations in the x Dimension');

% A Correlated Trajectory
x     = zeros(1,N);
c     = 0.80;           % degree of correlation; 0 to 1
step  = randn(1,N);
x(2) = randn();

for t=2:N
     x(t) = (c * x(t-1)) + ((1-c)*step(t));
end;

clf;
plot(xaxis, xcorr(x, 'coeff'));

xlabel('Lag');
ylabel('Correlation Coefficient');
title('Autocorrelation');
