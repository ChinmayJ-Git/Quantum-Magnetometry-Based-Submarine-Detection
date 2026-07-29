clear; clc; close all;
%rng(1);  

% 1. SENSOR ARRAY GEOMETRY  

armLen = 4;  % metres
sensorPos = [ 0,      0,      0;       % Sensor 1 (origin)
              armLen, 0,      0;       % Sensor 2
             -armLen, 0,      0;       % Sensor 4
              0,      armLen, 0;       % Sensor 5
              0,     -armLen, 0;       % Sensor 7
              0,      0,      armLen;  % Sensor 6
              0,      0,     -armLen]; % Sensor 3
nSensors = size(sensorPos, 1);


pairIdx = nchoosek(1:nSensors, 2);
nPairs  = size(pairIdx, 1);


% 2. AMBIENT GEOMAGNETIC FIELD DIRECTION  (Eq. 4: u = f(I, D))

I_deg = 60;   % inclination
D_deg = -9;   % declination
u = [cosd(I_deg)*cosd(D_deg);
     cosd(I_deg)*sind(D_deg);
     sind(I_deg)];


% 3. TRUE TARGET TRAJECTORY  

nSteps   = 21;
trueX    = linspace(-10, 10, nSteps);   % moves parallel to X-axis
trueY    = 10 * ones(1, nSteps);
trueZ    = -1 * ones(1, nSteps);
true_m   = 600;   % A.m^2   
true_th  = 60;    % degrees
true_phi = 15;    % degrees

noiseStd = 0.04;  % nT

% 4. PSO SETTINGS 

nParticles = 150;
maxIter    = 200;
c1 = 2; c2 = 2;
wMax = 0.9; wMin = 0.4;   % linearly decreasing inertia weight

lb = [-15, 5, -5, 0, -90];
ub = [15,15,5,90,90];


%5. MAIN TRACKING LOOP

estPose = zeros(nSteps, 5);   % [x y z theta phi] per step
estM    = zeros(nSteps, 1);
QI      = zeros(nSteps, 1);

for k = 1:nSteps
    targetPos = [trueX(k), trueY(k), trueZ(k)];

   
    dB = zeros(nSensors, 1);
    for s = 1:nSensors
        BA = dipoleField(sensorPos(s,:), targetPos, true_m, true_th, true_phi);
        dB(s) = u' * BA + noiseStd * randn();
    end

 
    dB_pairs = dB(pairIdx(:,1)) - dB(pairIdx(:,2));

    
    fitnessFcn = @(p) poseFitness(p, sensorPos, pairIdx, dB_pairs, u);

  
    if k==1
    initGuess=[trueX(1) trueY(1) trueZ(1) true_th true_phi];
else
    initGuess=estPose(k-1,:);
end

[bestP,~]=runPSO(fitnessFcn,lb,ub,nParticles,maxIter,c1,c2,wMax,wMin,initGuess);
    estPose(k, :) = bestP;

    
    f_ij = pairPredictions(bestP, sensorPos, pairIdx, u, 1);  
    valid = abs(f_ij) > 1e-6;                                  
    m_pairs = dB_pairs(valid) ./ f_ij(valid);

    m_pairs = m_pairs(isfinite(m_pairs));     
    m_pairs = m_pairs(abs(m_pairs) < 3000);    

    if isempty(m_pairs)
    estM(k) = NaN;
    QI(k) = 0;
    else
    estM(k) = median(m_pairs);
    QI(k) = qualityIndex(m_pairs);
end

    fprintf('Step %2d/%2d | true (x,y,z)=(%.1f,%.1f,%.1f)  est=(%.2f,%.2f,%.2f)  m_est=%.1f  QI=%.3f\n', ...
        k, nSteps, trueX(k), trueY(k), trueZ(k), bestP(1), bestP(2), bestP(3), estM(k), QI(k));
end

% 6. ERROR METRICS

rmseX = sqrt(mean((estPose(:,1)-trueX').^2));
rmseY = sqrt(mean((estPose(:,2)-trueY').^2));
rmseZ = sqrt(mean((estPose(:,3)-trueZ').^2));

rmsePos = sqrt(mean(sum((estPose(:,1:3) - [trueX' trueY' trueZ']).^2,2)));

fprintf('\n');
fprintf('RMSE X = %.4f m\n',rmseX);
fprintf('RMSE Y = %.4f m\n',rmseY);
fprintf('RMSE Z = %.4f m\n',rmseZ);
fprintf('Overall Position RMSE = %.4f m\n',rmsePos);



% 7. PLOTS

figure('Name', 'PSO Tracking Result', 'Color', 'w');

subplot(2,2,[1 3]);
plot3(trueX, trueY, trueZ, 'k-', 'LineWidth', 2); hold on;
plot3(estPose(:,1), estPose(:,2), estPose(:,3), 'ro--', 'MarkerFaceColor', 'r');
scatter3(sensorPos(:,1), sensorPos(:,2), sensorPos(:,3), 80, 'b', '^', 'filled');
grid on; axis equal; xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('True path', 'PSO estimate', 'Sensors', 'Location', 'best');
title('True vs. PSO-Estimated Target Trajectory');
view(35, 20);

subplot(2,2,2);
plot(1:nSteps, estM, 'b-o'); hold on;
yline(true_m, 'k--', 'True m');
xlabel('Time step'); ylabel('Estimated |m| (A m^2)');
title('Recovered Moment Magnitude'); grid on;

subplot(2,2,4);
plot(1:nSteps, QI, 'm-s');
xlabel('Time step'); ylabel('Quality Index');
ylim([0 1.05]); grid on;
title('Quality Index (accept/reject criterion)');

%  Estimated X vs. True X (parity plot) 
figure('Name', 'Estimated X vs True X', 'Color', 'w');
plot(trueX, trueX, 'k--', 'LineWidth', 1.2); hold on;   % y = x reference line
plot(trueX, estPose(:,1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; axis equal;
xlabel('True X (m)'); ylabel('Estimated X (m)');
title('PSO-Estimated  X vs. True X');
legend('Ideal (y = x)', 'PSO estimate', 'Location', 'best');



% LOCAL FUNCTIONS


function BA = dipoleField(sensorXYZ, targetXYZ, m, thetaDeg, phiDeg)

    r = sensorXYZ - targetXYZ;         r = -r;                        
    x = r(1); y = r(2); z = r(3);
    rnorm = norm(r);
    if rnorm < 1e-6
        rnorm = 1e-6;               
    end

    mu0 = 4*pi*1e-7;
    alpha = mu0*m / (4*pi*rnorm^3);
    beta  = 3*(x*cosd(thetaDeg)*cosd(phiDeg) + y*cosd(thetaDeg)*sind(phiDeg) + z*sind(thetaDeg)) / rnorm^2;

    Qvec = [cosd(thetaDeg)*cosd(phiDeg);
            cosd(thetaDeg)*sind(phiDeg);
            sind(thetaDeg)];

    BA = alpha * (beta*[x;y;z] - Qvec);  
    BA = BA * 1e9;                        
end

function f = pairPredictions(p, sensorPos, pairIdx, u, mUnit)

    x = p(1); y = p(2); z = p(3); th = p(4); phi = p(5);
    nSensors = size(sensorPos, 1);
    dBunit = zeros(nSensors, 1);
    for s = 1:nSensors
        BA = dipoleField(sensorPos(s,:), [x y z], mUnit, th, phi);
        dBunit(s) = u' * BA;
    end
    f = dBunit(pairIdx(:,1)) - dBunit(pairIdx(:,2));
end

function F = poseFitness(p, sensorPos, pairIdx, dB_pairs_meas, u)

    f_ij = pairPredictions(p, sensorPos, pairIdx, u, 1);
    denom = sum(f_ij.^2);
    if denom < 1e-9
        F = 1e12;   
        return;
    end
    m_hat = sum(dB_pairs_meas .* f_ij) / denom;   
    residual = dB_pairs_meas - m_hat * f_ij;
    F = sum(residual.^2);
end

function QI = qualityIndex(m_vals)

    n = numel(m_vals);
    if n < 2
        QI = 0;
        return;
    end
    s = 0;
    for pIdx = 1:n
        for qIdx = 1:n
            if pIdx ~= qIdx && abs(m_vals(qIdx)) > 1e-6
                s = s + ((m_vals(pIdx)/m_vals(qIdx)) - 1)^2;
            end
        end
    end
    QI = exp(-s * 10 / (n*(n-1)));   
end

function [gBest,gBestFit] = runPSO(fitnessFcn,lb,ub,nParticles,maxIter,c1,c2,wMax,wMin,initGuess)

    nDim = numel(lb);
    lb = lb(:)'; ub = ub(:)';

    % Initialization
    X = repmat(initGuess,nParticles,1) + 0.5*randn(nParticles,nDim);
    X = min(max(X,lb),ub);         % positions
    V = -(ub-lb) + 2*rand(nParticles, nDim).*(ub-lb);       % velocities
    V = V * 0.1;                                             % modest initial speed

    fit = zeros(nParticles, 1);
    for i = 1:nParticles
        fit(i) = fitnessFcn(X(i,:));
    end

    pBest = X;
    pBestFit = fit;
    [gBestFit, idx] = min(pBestFit);
    gBest = pBest(idx, :);

    %  Main PSO loop 
    for iter = 1:maxIter
        w = wMax - (wMax - wMin) * (iter / maxIter);   % linearly decreasing inertia

        for i = 1:nParticles
            r1 = rand(1, nDim);
            r2 = rand(1, nDim);
            V(i,:) = w*V(i,:) + c1*r1.*(pBest(i,:) - X(i,:)) + c2*r2.*(gBest - X(i,:));
            X(i,:) = X(i,:) + V(i,:);
            X(i,:) = min(max(X(i,:), lb), ub);   % keep particles in bounds

            fit(i) = fitnessFcn(X(i,:));
            if fit(i) < pBestFit(i)
                pBestFit(i) = fit(i);
                pBest(i,:)  = X(i,:);
            end
        end

        [minFit, idx] = min(pBestFit);
        if minFit < gBestFit
            gBestFit = minFit;
            gBest = pBest(idx, :);
        end
    end
end