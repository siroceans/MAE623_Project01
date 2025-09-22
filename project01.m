%% MAE 623 - CFD I: Project 01
% Jorge Luis Mares Zamora
% Due date: 09/23/2025

clear
clc
close all

%% Input parameters
scheme = 0; % Zero == Explicit; One == Implicit!!

alpha = 1; 
l = 1; 
h = 1; 
k = 1; 

To = 0;  % Initial temp. of the square
Tw = 0;  % West boundary - constant T bc
Tn = 0;  % North boundary - constant T bc
Tinf = 100;  % Freestream temperature

resolution_x = 10; 
resolution_y = 10; 

tfinal = 0.5; % set tfinal to 0 if you want a steady state solution!
Fo = 1/6; % make it small to guarantee stability!

t = 0; % Initial time
tol = 0.0001; % tolerance for steady state solution!

%% Creating the grid and initializing values

T = ones(resolution_x, resolution_y) * To; % creating grid where everything = To
Tnew = zeros(size(T)); 

dx = l / (resolution_x - 1); % square grid so in this case, dx == dy 
dy = l / (resolution_y - 1); 

% According to the handout Fo < 1/4 in 2D case
dt = Fo * dx^2 / alpha; % dt depending on our chosen value for Fo
Bi = h * dx / k; 

%% Explicit Time Scheme!

if scheme == 0
    while true
        % Calculating T at each interior node
        for m = 2:(resolution_x - 1) % Going through each row
            for n = 2:(resolution_y - 1) % and every column...
                Tnew(m,n) = Fo * (T(m+1, n) + T(m-1,n) + T(m,n+1) + T(m, n-1)) + (1 - 4 * Fo) * T(m, n); 
            end
        end

        % Applying boundary conditions (origin in matlab is top left!!)
        Tnew(:,1) = Tw; % west bc
        Tnew(1, :) = Tn; % north bc

        for n = 2:(resolution_y -1) % Insulated BC (T_south)

            m = 10; 
            Tnew(m, n) = Tnew(m-1, n); %T2 = T1
        end

        for m = 2:(resolution_x ) % Convective BC (T_east)
            n = 10; 
            Tnew(m, n) = (Bi * Tinf + Tnew(m, n-1)) / (1 + Bi); 
        end

        % Conditional to end the loop
        if tfinal == 0 % looking for steady state solution
            if max(max(abs(Tnew - T))) < tol
                break
            end
        elseif t >= tfinal % if we reached the max no of iterations
            break
        end

        t = t + dt; % update time
        T = Tnew;   % update T matrix for next iteration
    end
end

%disp(T)
plotFDResults(T, l, resolution_x, resolution_y, tfinal)

%% Implicit Time Scheme!

if scheme == 1
end

%% Plotting function
function plotFDResults(T, L, nx, ny, tfinal)
    % Plotting in 3D
    x = linspace(0, L, nx); 
    y = linspace(0, L, ny); 
    [X,Y] = meshgrid(x, y); 
    figure()
    surf(X, Y, flipud(T))
    title('2D Heat Conduction Equation.')
    xlabel('x direction')
    ylabel('y direction')
    zlabel('Temperature')

    % Plotting project requirements
    figure()
    subplot(2,1,1)
    yvalues = 0:(L/(ny-1)):L; 
    tvalues = T(:, nx/2)';
    plot(yvalues, tvalues)
    title('T vs. y for x = 0.5')
    ylabel('T values')
    xlabel('y values')
    if tfinal == 0
        legend('Steady State Solution', 'Location', 'best')
    else
        legend(['Solution at t = ', num2str(tfinal)])
    end

    subplot(2,1,2)
    xvalues = 0:(L/(nx-1)):L; 
    tvalues = T(ny/2, :); 
    plot(xvalues, tvalues)
    title('T vs. x for y = 0.5')
    ylabel('T values')
    xlabel('x values')
    if tfinal == 0
        legend('Steady State Solution', 'Location', 'best')
    else
        legend(['Solution at t = ', num2str(tfinal)])
    end
end