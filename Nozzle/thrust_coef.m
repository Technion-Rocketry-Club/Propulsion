clc; close all; clear
pa = 101325; % Atmospheric Pressure [Pa] (1 atm)
g = 1.215; % gamma
G = (2/(g + 1))^((g+1)/(2*(g-1))) * sqrt(g); % Big Gamma
t_b = 11; %[sec]

t = linspace(0,t_b,1001); % Time limit is middle value in seconds

d_t = 5; %[mm]
A_t = pi*(d_t/2000)^2; %[m^2]

% Ae/At ratio
% A = [1.25, 1.5, 1.75, 2, 2.25, 2.5]; % Can use to choose certain values
A = linspace(2, 20, 10); % linspace(X1, X2, N) generates N points between X1 and X2.
%A = linspace(1, 3, 15);

pres = [40, 50, 60]; % pc (chamber pressure) [bar]

for i = pres
    pci = i * 100000; % Initial Chamber Pressure [Pa]
    tau = 4.679; % [1/s]
    pc = pci .* exp(-t./tau); % Chamber Pressure as a function of time
    
    figure("Name", [num2str(i), ' bar'])
    hold on
    grid minor
       
    for j = 1:length(A)
        [mach,T,P,rho,area] = flowisentropic(g, A(j), 'sup'); % Returns pressure ratio: P
        pe = pc .* P; % Exit Pressure [Pa]
        cf = G .* sqrt((2.*g) ./ (g - 1) .* (1 - (pe./pc).^((g-1)./g))) + (pe./pc - pa./pc).*A(j); % Thrust Coeffecient 
    
        F = CalcThrust(cf, pc, A_t); %[N]
        I = CalcImpulse(F, t); %[N sec]

        plot(t, cf) 
        legendInfo{j} = ['A_e/A_t = ' num2str(A(j)) ', I = ' num2str(I,3) '[N sec]'];
    end

    xlabel('Time [sec]')
    ylabel('c_f')
    legend(legendInfo, 'Location', 'best')
    title("P_{ci} = " + pci/100000 + " [bar]")
    hold off
end

function F = CalcThrust(C_f, P_c, A_t)
    F = C_f .* P_c .* A_t; %[N]
end

function I = CalcImpulse(F, t)
    I = trapz(t, F); %[N sec]
end