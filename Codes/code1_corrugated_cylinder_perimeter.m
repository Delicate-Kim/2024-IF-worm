close all
clear all
clc

% Parameters for the corrugated cylinder

R = 0.25;      % Offset radius [mm]
A = 0.05*R;      % Amplitude of the corrugation, assume 5% of R
lambda = 0.5*R % Wavelength of the sine wave [mm], assume 50% of R
k = 2*pi/lambda; % Wave number
numPeriods = 13; %13; % L_a/lambda Number of periods (L_a = 1.625mm)
z_end = numPeriods * lambda; % End point of z for 10 periods

% Range of y_0 values from 0 to R + A
y0_vals = linspace(0, R + A, 1000);

% Pre-allocate for perimeter values (corrugated and normal cylinder)
perimeter_corrugated = zeros(size(y0_vals));
perimeter_normal = zeros(size(y0_vals));

% Define the function for dx/dz (derivative of x with respect to z) for corrugated cylinder
dx_dz = @(z, y0) (R + A .* sin(k.*z)) .* (A * k .* cos(k.*z)) ./ sqrt((R + A .* sin(k.*z)).^2 - y0^2);

% Loop over each y_0 and calculate the perimeter using numerical integration (corrugated cylinder)
for i = 1:length(y0_vals)
    y0 = y0_vals(i);
    
    % Corrugated cylinder perimeter
    if y0 <= R + A
        % Anonymous function for the integrand
        integrand = @(z) sqrt( (dx_dz(z, y0)).^2 + 1 );

        % Set the limits for z (from 0 to z_end)
        z1 = 0;
        z2 = z_end;

        % Perform numerical integration using integral function
        %perimeter_corrugated(i) = integral(integrand, z1, z2);
        perimeter_corrugated(i) = 4 * sqrt(R^2 - y0^2) + 2*integral(integrand, z1, z2);
    else
        perimeter_corrugated(i) = NaN; % For values of y_0 beyond the physical boundary
    end
    
    % Normal cylinder perimeter
    if y0 <= R
        %perimeter_normal(i) = 2 * pi * sqrt(R^2 - y0^2);
        perimeter_normal(i) = 4 * sqrt(R^2 - y0^2) + 2*z_end;
    else
        perimeter_normal(i) = NaN; % For values of y_0 beyond the radius of the normal cylinder
    end
end

% Plot the perimeter as a function of y_0 for both cylinders
figure;
plot(y0_vals, perimeter_corrugated, 'Color', [0 0.5 0.8],  'LineWidth', 3); hold on;
plot(y0_vals, perimeter_normal, 'Color', [0.9 0.3 0], 'LineWidth', 3);
xlabel('Intersecetion Height (mm)', 'FontSize', 12);
ylabel('Intersection Perimeter (mm)', 'FontSize', 12);
%title('Intersection perimeter vs height', 'FontSize', 14);
legend('Corrugated', 'Smooth', 'Location', 'NorthWest');

set(gca, 'FontSize', 20)

set(gcf,'position',[0,0,700,500])
ylim([2 8])
%xline(0.21, 'k.');
%xline(0.22, 'k.');
%xline(0.23, 'k.');
%xline(0.24, 'k.');
%xline(0.25, 'k.');
grid on;
