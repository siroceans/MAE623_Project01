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

tfinal = 0; % set tfinal to 0 if you want a steady state solution!
Fo = 1/6; % make it small to guarantee stability!

tol = 0.0001; % tolerance for steady state solution!


%% Functioning work
% Uncomment this section if you want to modify the code using the input parameters for resolution x and resolution y
% returns the T matrix of the solution

%explicitFD(To, resolution_x, resolution_y, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k)
%implicitFD(To, resolution_x, resolution_y, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k) 

%% Plotting Results for Report

T10E = explicitFD(To, 10, 10, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k); 
T20E = explicitFD(To, 20, 20, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k); 
T40E = explicitFD(To, 40, 40, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k); 



T10I = implicitFD(To, 10, 10, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k);
T20I = implicitFD(To, 20, 20, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k);
T40I = implicitFD(To, 40, 40, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k);

% Plotting Explicit y
y1 = 0:(l/9):l; 
T1 = T10E(:, 5); 
y2 = 0:(l/19):l; 
T2 = T20E(:,10); 
y3 = 0:(l/39):l;
T3 = T40E(:, 20); 
figure() 
plot(y1,T1, y2,T2, y3,T3)
ylabel('T values')
xlabel('y values')
legend('10x10 grid', '20x20 grid', '40x40 grid', 'Location', 'Best')
if tfinal > 0
    title('T vs. y for x = 0.5 at t = 0.05 (Explicit)')
else 
    title('T vs. y for x = 0.5 at steady state (Explicit)')
end

% Plotting Explicit x
x1 = y1; 
x2 = y2; 
x3 = y3; 
T1 = T10E(5, :); 
T2 = T20E(10, :); 
T3 = T40E(20, :); 
figure()
plot(x1, T1, x2, T2, x3, T3)
xlabel('x values')
ylabel('T values')
legend('10x10 grid', '20x20 grid', '40x40 grid', 'Location', 'Best')
if tfinal > 0
    title('T vs. x for y = 0.5 at t = 0.05 (Explicit)')
else 
    title('T vs. x for y = 0.5 at steady state (Explicit)')
end

% Plotting Implicit y
T1 = T10I(:, 5); 
T2 = T20I(:, 10); 
T3 = T40I(:, 20);
figure()
plot(y1,T1, y2,T2, y3,T3)
ylabel('T values')
xlabel('y values')
legend('10x10 grid', '20x20 grid', '40x40 grid', 'Location', 'Best')
if tfinal > 0
    title('T vs. y for x = 0.5 at t = 0.05 (Implicit)')
else 
    title('T vs. y for x = 0.5 at steady state (Implicit)')
end

% Plotting Implicit x
T1 = T10I(5, :); 
T2 = T20I(10, :); 
T3 = T40I(20, :); 
figure()
plot(x1, T1, x2, T2, x3, T3)
xlabel('x values')
ylabel('T values')
legend('10x10 grid', '20x20 grid', '40x40 grid', 'Location', 'Best')
if tfinal > 0
    title('T vs. x for y = 0.5 at t = 0.05 (Implicit)')
else 
    title('T vs. x for y = 0.5 at steady state (Implicit)')
end



%% Explicit Time Scheme Function!

function T = explicitFD(To, resolution_x, resolution_y, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k)

    % Creating Grid and Initializing Values
    dx = l / (resolution_x - 1); % square grid so in this case, dx == dy 

    % According to the handout Fo < 1/4 in 2D case
    dt = Fo * dx^2 / alpha; % dt depending on our chosen value for Fo
    Bi = h * dx / k; 

    t = 0; % Initial time;
    T = ones(resolution_x, resolution_y) * To; % creating grid where everything = To
    Tnew = zeros(size(T)); 
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

            m = resolution_x; 
            Tnew(m, n) = Tnew(m-1, n); %T2 = T1
        end

        for m = 2:(resolution_x ) % Convective BC (T_east)
            n = resolution_y; 
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
        T = Tnew;  % update T matrix for next iteration
    end
    %disp(t)
%    plotFDResults(T, l, resolution_x, resolution_y, tfinal)
end

%% Implicit Time Scheme Function!

function T = implicitFD(To, resolution_x, resolution_y, Fo, tfinal, Tw, Tn, Tinf, tol, l, alpha, h, k)

    % Creating Grid and Initializing Values
    dx = l / (resolution_x - 1); % square grid so in this case, dx == dy 

    % According to the handout Fo < 1/4 in 2D case
    dt = Fo * dx^2 / alpha; % dt depending on our chosen value for Fo
    Bi = h * dx / k; 

    t = 0; % Initial time
    % Creating column matrices 
    T = ones(resolution_x * resolution_y, 1) * To; % column matrix for initial T
    Tnew = zeros(resolution_x * resolution_y, 1); % column matrix for T!
    N = resolution_x * resolution_y; 

    % Creating matrix A
    dp = 1 + 4 * Fo; % Point coefficient values
    dn = - Fo; % Neighbor coefficient valuesa

    % Array of Diagonal Values
    diagonalArray = zeros(N, 5);
    diagonalArray(:, 1) = dn; 
    diagonalArray(:, 2) = dn;
    diagonalArray(:, 3) = dp; 
    diagonalArray(:, 4) = dn; 
    diagonalArray(:, 5) = dn;

    diagonalPositions = [-resolution_y, -1, 0, 1, resolution_y]; % Locations of the diagonals!
    A = spdiags(diagonalArray, diagonalPositions, N, N); % Creating the sparce matrix

    % Modifying rows for boundary conditions!
    left = 1 : resolution_y; 
    right = resolution_y^2 - resolution_y  + 1: resolution_y^2; 
    bottom = (0:resolution_x - 1) * resolution_x + 1;  
    top = (1: resolution_x) * resolution_x; 

    for p = 1:length(T)
        if ismember(p, left)
            % left BC
            A(p, :) = 0;
            A(p, p) = 1; 
        elseif ismember(p, top)
            % top BC
            A(p, :) = 0;
            A(p, p) = 1; 
        elseif ismember(p, right)
            % right BC
            A(p, :) = 0;
            A(p, p) = (1 + Bi); 
            A(p, p - resolution_y) = -1; 
        elseif ismember(p, bottom)
            % bottom BC
            A(p, :) = 0;
            A(p, p) = -1; 
            A(p,p + 1) = 1; 
        end
    end

    while true % time marching loop
        % Creating C vector
        C = T; 
        for p = 1:length(C) % Applying BC's
            if ismember(p, left)
                C(p) = Tw; 
            elseif ismember(p, top)
                C(p) = Tn; 
            elseif ismember(p, right)
                C(p) = Bi * Tinf; 
            elseif ismember(p, bottom)
                C(p) = 0; 
            end
        end

        Tnew = A \ C; % Solving equation for Tnew

        % Conditional to end the loop
        if tfinal == 0 % looking for steady state sol'n
            if max(abs(Tnew - T)) < tol
                break 
            end 
        elseif t >= tfinal  % Reached ending t value
            break 
        end

        t = t + dt; 
        T = Tnew; % Updating T vector for new iteration
    end
    T = reshape(T, [resolution_x, resolution_y]); 
    T = flipud(T); 
%    plotFDResults(T, l, resolution_x, resolution_y, tfinal)
    %disp(t)
end

%% Plotting function
function plotFDResults(T, L, nx, ny, tfinal)
    % Flipping Matrix!
    T = flipud(T); 

    % Plotting in 3D
    x = linspace(0, L, nx); 
    y = linspace(0, L, ny); 
    [X,Y] = meshgrid(x, y); 
    figure()
    surf(X, Y, T)
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