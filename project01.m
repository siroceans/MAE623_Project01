%% MAE 623 - CFD I: Project 01
% Jorge Luis Mares Zamora
% Due date: 09/23/2025

clear
clc

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

t = 0; % Initial timea
tfinal = 0.5; % set tfinal to 0 if you want a steady state solution!

%% Creating the grid and initializing values

T = ones(resolution_x, resolution_y) * To; % creating grid where everything = To
dx = l / (resolution_x - 1); % square grid so in this case, dx == dy 
dy = l / (resolution_y - 1); 

% According to the handout Fo < 1/4 in 2D case
% I will choose Fo = 1/2 to guarantee stability!
Fo = 1/2; 
dt = Fo * dx^2 / alpha; % dt depending on our chosen value for Fo
Bi = h * dx / k; 

%% Explicit Time Scheme!

Tnew = zeros(size(T)); 

if scheme == 0
    while true
        % Calculating T at each interior node
        for m = 2:(resolution_x - 1) % Going through each row
            for n = 2:(resolution_y - 1) % and every column...
                Tnew(m,n) = Fo * (T(m+1, n) + T(m-1,n) + T(m,n+1) + T(m, n-1)) + (1 - 4 * Fo) * T(m, n); 
%                fprintf('m is %d, n is %d \n', m, n)
            end
        end

        % Applying boundary conditions (origin in matlab is top left!!)
        Tnew(:,1) = Tw; % west bc
        Tnew(1, :) = Tn; % north bc

        % Insulated BC (T_south)
        for n = 2:(resolution_y -1)
            m = 10; 
            Tnew(m, n) = Fo * (2 * T(m-1, n) + T(m, n+1) + T(m, n-1)) + (1 - 4 * Fo) * T(m, n); 
        end

        % Convective BC (T_east)
        for m = 2:(resolution_x )
            n = 10; 
            Tnew(m, n) = (Bi * Tinf + Tnew(m, n-1)) / (1 * Bi); 
        end

        disp(Tnew)
        break
    end
end