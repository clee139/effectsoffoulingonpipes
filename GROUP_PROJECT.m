% PROJECT #2
% Name: Jaylin Trice
% ME 2543--Simulation Methods
% SPRING 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; clear all% Clear the variable list and the command window

% PROBLEM #: 1

% DEFINE THE GLOBAL VARIABLES to BE USED IN CB
global e1 D rho L mu;
% Declare values as global variables to use them inside the function
global L1 L2 L3 D1 D2 D3 A1 A2 A3 Q_Total;
global j1 j2 j3;

% DEFINE THE GIVEN VALUES THAT DON'T CHANGE
rho = 62.4179721; % density of water lbm/ft^3; 32.174
mu = 2.05e-5; % units lbf*s/ft^2
e1 = 0.00085; % surface roughness; unitless
Q_Total = 1.67101; % ft^3/s

% DEFINE THE PIPE DIMENSIONS
L1 = 10; % length; feet
L2 = L1;
L3 = L1;

D1 = 2.0/12; % Diameter; feet
D2 = D1;
D3 = D1;

A1 = pi()*(1/4)*(D1)^2; % ft^2
A2 = A1;
A3 = A1;

j1 = e1/D1;
j2 = e1/D2;
j3 = e1/D3;

% initial guess of f; NOTE: I SOLVED FOR f THRU BRUTE FORCE FIRST
f0 = 0.0058; % estimating from EQ(9)

% solve for Q1, Q2, Q3
handle = @myH;

% chose random initial guess for Q1-Q3; chose 5 for pipe
% because there shouldn't be a huge pressure drop
X0 = [0.2 0.2 0.2 30 0.005 0.005 0.005]; % [x(1) x(2) x(3) x(4)]
answer = fsolve(handle, X0);
disp("The flow rate for pipe 1 is : " + answer(1)*448.831 + " gpm")
disp("The flow rate for pipe 2 is : " + answer(2)*448.831 + " gpm")
disp("The flow rate for pipe 3 is : " + answer(3)*448.831 + " gpm")
disp("The pressure drop is : " + answer(4)*0.00691964 + " lbf/in^2")

% factor 1: 0.00021583993174700825
% factor 2: 0.006944444

%% PROBLEM 2 %%

% main issue:
% figure out why I don't get same values as myf(x)
% figure out how to calculate Reynolds number with given info
% figure out if the laminar or turbulent flow

%the individual pipe properties from table 1
L1 = 20; % feet
L2 = 10; % feet
L3 = 30; % feet
D1 = 2.0/12; % feet
D2 = 2.5/12; % feet
D3 = 1.5/12; % feet
A1 = pi()*(1/4)*(D1)^2; % ft^2
A2 = pi()*(1/4)*(D2)^2; % ft^2
A3 = pi()*(1/4)*(D3)^2; % ft^2
j1 = e1/D1;
j2 = e1/D2;
j3 = e1/D3;

% answer if turbulent flow
test_2 = fsolve(@myH, [0.3 0.3 0.3 100 0.005 0.005 0.005]);
% answer if laminar flow
test_3 = fsolve(@myG, [0.3 0.3 0.3 100 0.005 0.005 0.005]);
pressure_drop = test_2(4)*0.00691964;
pressure_drop2 = test_3(4)*0.00691964;

% the values for turbelent flow look more reasonable so use them
disp("The flow rate for pipe 1 is : " + test_2(1)*448.831 + " gpm")
disp("The flow rate for pipe 2 is : " + test_2(2)*448.831 + " gpm")
disp("The flow rate for pipe 3 is : " + test_2(3)*448.831 + " gpm")
disp("The pressure drop is : " + pressure_drop + " lbf/in^2")

% rho, mu and e values are all the same
% since pipe dimensions are variable we must calculate Re for each pipe

%% PROBLEM 3

% use a for loop to go through all the Q values; incremenets of 100
Q_Values = 100:100:1500;
Q_imperial = Q_Values/448.831;

ans_laminar = [];
ans_turbulent = [];
X0 = [0.3 0.3 0.3 100 0.005 0.005 0.005];
pressureDrop = [];

for i = 1:length(Q_Values)
    Q_Total = Q_imperial(i);
    ans_laminar = [ans_laminar; fsolve(@myG, X0)];
    ans_turbulent = [ans_turbulent; fsolve(@myH, X0)];
end

% use for loop to get all the pressure drops from each row in ans_turbulent
for i = 1:length(ans_turbulent)
    p = ans_turbulent(i,4)*0.00691064;
    pressureDrop = [pressureDrop p];
end

%% PROBLEM 4

% part A: fouling so e/D increases by 25%
j1 = (e1/D1)*1.25;
j2 = (e1/D2)*1.25;
j3 = (e1/D2)*1.25;

answer_laminar = [];
answer_turbulent = [];
pDrop = [];

for i = 1:length(Q_Values)
    Q_Total = Q_imperial(i);
    answer_laminar = [answer_laminar; fsolve(@myG, X0)];
    answer_turbulent = [answer_turbulent; fsolve(@myH, X0)];
end

% use for loop to get all the pressure drops from each row in ans_turbulent
for i = 1:length(answer_turbulent)
    p = answer_turbulent(i,4)*0.00691964;
    pDrop = [pDrop p];
end

% part B: fouling where e/D increases by 35%
j1 = (e1/D1)*1.35;
j2 = (e1/D2)*1.35;
j3 = (e1/D2)*1.35;

answer_laminar2 = [];
answer_turbulent2 = [];
pDrop2 = [];

for i = 1:length(Q_Values)
    Q_Total = Q_imperial(i);
    answer_laminar2 = [answer_laminar2; fsolve(@myG, X0)];
    answer_turbulent2 = [answer_turbulent2; fsolve(@myH, X0)];
end

% use for loop to get all the pressure drops from each row in ans_turbulent
for i = 1:length(answer_turbulent2)
    p = answer_turbulent2(i,4)*0.00691964;
    pDrop2 = [pDrop2 p];
end

%% Plot answers from Parts 3 and 4

% plot(Q_Values,pressureDrop,'Color','r')
% hold on
% plot(Q_Values,pDrop,'Color','b')
% plot(Q_Values,pDrop2,'Color','g')
% title("Total Flow Rate vs. Pressure Drop")
% legend("No fouling", "25% fouling", "35% fouling")
% xlabel("Total Flow Rate (gpm)")
% ylabel("Pressure Drop (psi)")

%% FUNCTIONS %%

% function for variable L and D values
% x(1)-(3) are the Q values; x(4) is the dP value; x(5)-(7) are f values
function [F] = myG(x) % assumes laminar flow
    global D1 D2 D3 L1 L2 L3 rho mu Q_Total
    A1 = pi()*D1/4;
    A2 = pi()*D2/4;
    A3 = pi()*D3/4;

    F = [ x(5) *(L1/D1)*(1/2)*((x(1)/A1)^2)-(x(4)/rho);
          x(6) *(L2/D2)*(1/2)*((x(2)/A2)^2)-(x(4)/rho);
          x(7) *(L3/D3)*(1/2)*((x(3)/A3)^2)-(x(4)/rho);
          x(1) + x(2) + x(3) - Q_Total;
          x(5) - 64/(rho*x(1)*D1/(mu*A1));
          x(6) - 64/(rho*x(2)*D2/(mu*A2));
          x(7) - 64/(rho*x(3)*D3/(mu*A3))];
end

function [F] = myH(x) % assumes turbulent flow
    global D1 D2 D3 L1 L2 L3 rho mu e1 Q_Total j1 j2 j3
    A1 = pi()*(D1)/4;
    A2 = pi()*(D2)/4;
    A3 = pi()*(D3)/4;
    
    F = [ x(5)*(L1/D1)*(1/2)*((x(1)/A1)^2)-(x(4)/rho);
          x(6)*(L2/D2)*(1/2)*((x(2)/A2)^2)-(x(4)/rho);
          x(7)*(L3/D3)*(1/2)*((x(3)/A3)^2)-(x(4)/rho);
          x(1) + x(2) + x(3) - Q_Total;
          (1/sqrt(x(5))) + 2 * log(((j1)/3.7)+(2.51/((rho*x(1)*D2/(mu*A1))*sqrt(x(5)))));
          (1/sqrt(x(6))) + 2 * log(((j2)/3.7)+(2.51/((rho*x(2)*D2/(mu*A2))*sqrt(x(6)))));
          (1/sqrt(x(7))) + 2 * log(((j3)/3.7)+(2.51/((rho*x(3)*D3/(mu*A3))*sqrt(x(7)))))];
end