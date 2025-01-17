%% This code plot dimensionless net force of 2D cylinder partially immersed in water for different Bond number and Specific gravity. The code also plot theoretical prediction of meniscus angle and fit a polynomial function.
%% Writer: Soohwan Kim, Georgia Tech ME, ksh35000@gatech.edu

%% calculate phi under equilibrium, F_net = 0 using fsolve, initial estimation: z = l_c * tan(alpha) -> better estimation: z/R = sqrt((2/Bo)*(1+cos(theta+phi))) = sqrt((2/Bo)*(1-cos(alpha)))

clear all;
close all;
clc

%% test different combination - worm values
theta = 90*pi/180; % [rad]
%alpha = 30*pi/180;
%phi = pi - theta + alpha;
sigma = 72; %72; %sigma = 69.3; % [dyne/cm]  
rho_f = 1; % [g/cm^3]
g = 980; % [cm/s^2]
l_c = sqrt(sigma/rho_f/g) % [cm]
R = 0.25;  %0.023; %R = 0.2633 pc; %R = 0.26795; epoxy % [cm]
%%

Bo = rho_f*g*R^2/sigma; % Bo = 0.85


figure(1)
set(gca,'FontSize',18)
xlabel('Water curvature angle (degree)');
ylabel('Dimensionless F_{net}');
xlim([0 360])
grid on;
hold on;

yline(0, 'k--', 'LineWidth', 1)

% R = 0.025
SG = 1.0;
%fplot (@(a) 2*sin(a) + Bo*(2*l_c*tan(a)/R*sin(pi+a-theta) + (pi+a-theta) - 0.5*sin(2*pi+2*a-2*theta) - SG*pi),[0 2*pi],'k');
p = fplot (@(a) 2*sind(a) + Bo*(2*sqrt((2/Bo)*(1-cosd(a)))*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'b', 'LineWidth', 2); % a in degree
%p = fplot (@(a) 2*sind(a) + Bo*(2*l_c*tand(a)/R*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2); %old
text2line(p, 0.454, 0,'SG = 1.0')
SG = 1.5;
p = fplot (@(a) 2*sind(a) + Bo*(2*sqrt((2/Bo)*(1-cosd(a)))*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2);
%p = fplot (@(a) 2*sind(a) + Bo*(2*l_c*tand(a)/R*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2); %old
text2line(p, 0.5, 0,'SG = 1.5')
SG = 2;
p = fplot (@(a) 2*sind(a) + Bo*(2*sqrt((2/Bo)*(1-cosd(a)))*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2);
%p = fplot (@(a) 2*sind(a) + Bo*(2*l_c*tand(a)/R*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2); %old
text2line(p, 0.5, 0,'SG = 2.0')
SG = 2.5;
p = fplot (@(a) 2*sind(a) + Bo*(2*sqrt((2/Bo)*(1-cosd(a)))*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2);
%p = fplot (@(a) 2*sind(a) + Bo*(2*l_c*tand(a)/R*sind(pi*180/pi+a-theta*180/pi) + (pi+a*pi/180-theta) - 0.5*sind(2*pi*180/pi+2*a-2*theta*180/pi) - SG*pi),[0 360],'k', 'LineWidth', 2); %old
text2line(p, 0.5, 0,'SG = 2.5')

xticks([0 90 180 270 360])

set(gcf,'position',[500,0,600,500])

%%
SG = 1:0.01:2;% R = 0.25
%SG = 1:0.5:100; % R = 0.025

syms a
alpha_arr1 = [];
alpha_arr2 = [];
alpha_arr3 = [];
alpha_arr4 = [];
alpha_arr5 = [];

for i=1:length(SG)
    equation = @(a) 2*sin(a) + Bo*(2*sqrt((2/Bo)*(1-cos(a)))*sin(pi+a-theta) + (pi+a-theta) - 0.5*sin(2*pi+2*a-2*theta) - SG(i)*pi) - 0; % a in radian
    %sol = solve(2*sin(a) + Bo*(2*l_c*tan(a)/R*sin(pi+a-theta) + (pi+a-theta) - 0.5*sin(2*pi+2*a-2*theta) - SG(i)*pi)==0,a);
    [sol1, fval1, exitflag1] = fsolve(equation, 0);
    [sol2, fval2, exitflag2] = fsolve(equation, 100*pi/180);
    [sol3, fval3, exitflag3] = fsolve(equation, 270*pi/180);
    if exitflag1 > 0
        alpha_arr1 = [alpha_arr1 sol1*(180/pi)];
    else
        alpha_arr1 = [alpha_arr1 -360]; % give -360 value for non-existing solution case
    end

    if exitflag2 > 0
        alpha_arr2 = [alpha_arr2 sol2*(180/pi)];
    else
        alpha_arr2 = [alpha_arr2 -360];
    end

    if exitflag3 > 0
        alpha_arr3 = [alpha_arr3 sol3*(180/pi)];
    else
        alpha_arr3 = [alpha_arr3 -360];
    end

end

figure(2)
hold on;
size = 22.5;
%alpha_arr1(137) = -360;
%alpha_arr1(138) = -360;
%alpha_arr1(144) = -360
plot(SG, alpha_arr1, '.','MarkerFaceColor', '#2B96D4','MarkerSize', size)
%alpha_arr2(137) = -360;
%alpha_arr2(138) = -360;
%alpha_arr2(144) = -360
plot(SG, alpha_arr2, '.','MarkerFaceColor', '#E64B00', 'MarkerSize', size)
plot(SG, alpha_arr3, '.','MarkerFaceColor', '#4CAF50', 'MarkerSize', size)

xlabel('SG');
ylabel('Water curvature angle (degree)');
ylim([0 360])
%yticks([0 90 180 270 360])
grid on;
set(gca,'FontSize',18)

set(gcf,'position',[1000,0,600,500])

%sol = solve(2*sin(a) + Bo*(2*l_c*tan(a)/R*sin(pi+a-theta) + (pi+a-theta) - 0.5*sin(2*pi+2*a-2*theta) - SG*pi)==0,a);

%Z = 2*sin(alpha) + X.*(2*sin(phi)*z_c/R + phi - 0.5*sin(2*phi) - Y.*pi);

%SG_harry=[1.056222979 1.416433505 1.158265621 1.209286942 1.260308263 1.311329584...
%1.362350905 1.413372226 1.464393547 1.489904208 1.49500634 1.500108472 1.515414868...
%1.566436189 1.61745751 1.642968171 1.648070303 1.653172435 1.658274567 1.663376699...
%1.664397126 2.331756004 3.60728903];

%alpha_harry=[14.17831672 29.50268848 17.53047202 19.63426428 21.83891572 24.16719474...
%26.65044662 29.33395962 32.28780936 33.90017978 34.23601777 34.57679399 35.63124419...
%39.60104343 44.83180245 48.65161766 49.66393955 50.85097316 52.35327193 54.80273823...
%57.6083509 50.85097316 50.85097316];

%plot(SG_harry, alpha_harry, 'r.','MarkerSize', 10)


yline(90, 'k--', 'LineWidth', 1)
xlim([1 2])
%xlim([1 100])
%% close up
figure(4)
hold on;
size = 25;
plot(SG, alpha_arr1, 'k.','MarkerSize', size)

xlabel('SG');
ylabel('Water curvature angle (degree)');
xlim([1 2])
ylim([0 90])
yticks([0 15 30 45 60 75 90])
grid on;
set(gca,'FontSize',18)

%% polynomial fitting

% Extracting SG and alpha values from figure(4)
SG_closeup = SG;
alpha_closeup = alpha_arr1;

% Remove any -360 values indicating non-existing solutions
valid_indices = alpha_closeup ~= -360;
SG_valid = SG_closeup(valid_indices);
alpha_valid = alpha_closeup(valid_indices);

% Perform polynomial fitting
% Here, a 2nd-degree polynomial is chosen, but you can adjust the degree as needed
degree = 10; % Change degree for a higher/lower order polynomial
[p, S] = polyfit(SG_valid, alpha_valid, degree);

% Generate fitted values for plotting
SG_fitted = linspace(min(SG_valid), max(SG_valid), 100);
alpha_fitted = polyval(p, SG_fitted);

% Plot the original data and the polynomial fit
figure(5)
hold on;
plot(SG_valid, alpha_valid, 'ko', 'MarkerSize', 5) % Original data points
plot(SG_fitted, alpha_fitted, 'r-', 'LineWidth', 2) % Polynomial fit

xlabel('SG');
ylabel('Water curvature angle (degree)');
title('Polynomial Fit for SG vs Water Curvature Angle');
grid on;
set(gca, 'FontSize', 18);
legend('Original Data', ['Polynomial Fit (Degree = ' num2str(degree) ')'], 'Location', 'Best');

% Display polynomial coefficients
disp('Polynomial coefficients:');
disp(p);

% Display the polynomial equation
equation = ['y = ' num2str(p(1)) ' * x^' num2str(degree)];
for i = 2:degree
    equation = [equation ' + ' num2str(p(i)) ' * x^' num2str(degree-i+1)];
end
equation = [equation ' + ' num2str(p(end))];
disp('Polynomial Equation:');
disp(equation);

%%

R = 0.2825; % [cm]
L = 8; % [cm]
V = R^2*pi*L; %[cm^3]

m_straw = 208.2e-3; % [g] 
m_added = 1317e-3; % [g]
m = 3301.
%% tan(alpha) vs sqrt(2*(1-cos(alpha)))
%figure(3)
%hold on
%fplot(@(a) tan(a) , 'k')
%fplot(@(a) sqrt(2*(1-cos(a))))
%xlim([0 2*pi])

function text2line(h,ksi,z,T)
% Inserts text T in/near line with handle h
%  ksi - relative distance from the beginning of curve,
%  z - shift along normal to curve
%
set(gcf, 'CurrentObject', h)
x=h.XData;
y=h.YData;
i = round(ksi*numel(x));
% Get the local slope
dy=y(i+1)-y(i-1);
dx=x(i+1)-x(i-1);
d = dy/dx;
X = diff(get(gca, 'xlim'));
Y = diff(get(gca, 'ylim'));
p = pbaspect;
a = atan(d*p(2)*X/p(1)/Y)*180/pi;
% Display the text
switch z==0
    case 1
        text(x(i), y(i), T,'HorizontalAlignment','center', 'BackgroundColor', 'w', 'rotation', a, 'FontSize', 10);
    case 0
        ez=[dy,-dx]/norm([dy,-dx]); % unit normal vector
        text(x(i)+z*ez(1), y(i)+z*ez(2), T, 'HorizontalAlignment','center', 'rotation', a, 'FontSize', 10);
end

end