%% plot a radial gradient, zoom in on part of it

rmax = 100;
r = linspace(0, 1, rmax);                               % Define Radius & Radius Gradient Vector
a = linspace(0, 2*pi, 150);                             % Angles (Radians)
[R,A] = ndgrid(r, a);                                   % Create Grid
Z = R;   
Z2 = A;                 
[X,Y,Z] = pol2cart(A,R,Z);   
[X,Y,Z2] = pol2cart(A,R,Z2);                              % Convert To Cartesian

figure(1)
surf(X, Y, Z)
view(0, 90)
axis('equal')
colormap(cmocean('tempo'))
colorbar('Ticks', [])
shading('interp')
axis([-0.7 0 -0.7 0])
set(gca, 'Visible','off')

%% threshold it

indices = Z>0.5;
thresholdZ = Z;
thresholdZ(indices) = 1;
thresholdZ(~indices) = 0;

figure(2)
surf(X, Y, thresholdZ)
view(0, 90)
axis('equal')
colormap(cmocean('tempo'))
shading('interp')
axis([-0.7 0 -0.7 0])
set(gca, 'Visible','off')

%% band it

indices = Z > 0.45 & Z < 0.5;
thresholdZ = Z;
thresholdZ(indices) = 0;
thresholdZ(~indices) = 1;

figure(3)
surf(X, Y, thresholdZ)
view(0, 90)
axis('equal')
colormap(cmocean('gray'))
shading('interp')
axis([-0.7 0 -0.7 0])
set(gca, 'Visible','off')

%% band but resistance colored

indices = Z >= 0.45 & Z <= 0.5;
thresholdZ2 = Z2;
thresholdZ2(~indices) = nan; %1.25*pi;

figure(4)
surf(X, Y, thresholdZ2)
view(0, 90)
axis('equal')
colormap(cmocean('-balance'))
caxis([pi pi*1.5])
shading('interp')
axis([-0.7 0 -0.7 0])
set(gca, 'Visible','off')