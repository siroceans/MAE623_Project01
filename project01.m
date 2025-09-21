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

t = 0; % Initial time

%% Creating the grid

T = ones(resolution_x, resolution_y) * To; % creating grid where everything = To
dx = l / (resolution_x - 1); % square grid so in this case, dx == dy 
dy = l / (resolution_y - 1); 

% According to the handout Fo < 1/4 in 2D case
% I will choose Fo = 1/2 to guarantee stability!
Fo = 1/2; 
dt = Fo * dx^2 / alpha; % dt depending on our chosen value for Fo

