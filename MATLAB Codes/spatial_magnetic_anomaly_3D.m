% Simulates the spatial magnetic anomaly produced by a magnetic dipole
% using the full magnetic dipole equation.

clear;
clc;
close all;

% Physical constants

mu0 = 4*pi*1e-7;             


m = [0; 0; 1e6];              % Dipole aligned along +z

% Observation plane

x = linspace(-800,800,250);
y = linspace(-800,800,250);

[X,Y] = meshgrid(x,y);

z0 = 100;                     % Aircraft altitude (m)

% Preallocate

Bx = zeros(size(X));
By = zeros(size(X));
Bz = zeros(size(X));

% Compute magnetic field

for i = 1:size(X,1)

    for j = 1:size(X,2)

        r = [X(i,j); Y(i,j); z0];

        rNorm = norm(r);

        B = (mu0/(4*pi)) * ...
            ( (3*r*(dot(m,r))/rNorm^5) - (m/rNorm^3) );

        Bx(i,j) = B(1);
        By(i,j) = B(2);
        Bz(i,j) = B(3);

    end

end

% Magnitude of magnetic anomaly

Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);


Bmag = Bmag*1e9;

% Plot

figure('Color','w');

surf(X,Y,Bmag,...
    'EdgeColor','none');

shading interp

colormap(parula)

cb = colorbar;
cb.Label.String = 'Magnetic Anomaly, B_A (nT)';
cb.FontSize = 11;

xlabel('x / m','FontSize',12)
ylabel('y / m','FontSize',12)
zlabel('Magnetic Anomaly, B_A (nT)','FontSize',12)

title('Spatial Distribution of Magnetic Anomaly','FontSize',14)

view(45,30)

grid on
box on

set(gca,...
    'FontSize',11,...
    'LineWidth',1.2);

axis tight