% SimulateParticles Function

function particle = SimulateParticles(N, particleCount, tau, k, dimensions)

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

end