%% Question01
clc
clear all
warning off

G1 = tf([1,1], [1,2]);
G2 = tf([1], [500, 0, 0]);
G = series(G1, G2);
unity = tf([1],[1]);
system = feedback(G, unity, -1);

% Trasnfer Function
system

% Poles
poles = pole(system)

% Pole-Zero Plot
figure
pzmap(system)
xlim([-3, 1.5])
ylim([-0.08, 0.08])

% Is Stable?
isStable = isstable(system)

% Impulse Function
figure
impulse(system)

% Bode Plot
figure
bode(system)

% Output
figure
t = 0:0.01:100000;
u = exp(-t./10000).*cos(pi/500*t);
lsim(system, u, t)

%% Question02
%% 2.1
wn = 10;
table = [5, 1, 0.5, 0, -0.5];
for mu = table
    system = tf(wn^2,[1,2*mu*wn,wn^2]);
    bode(system)
    hold on
end
legend('\mu = 5','\mu = 1','\mu = 0.5','\mu = 0','\mu = -0.5');

figure
for mu = table
    if (mu == -0.5)
        xlim([0,2.5]);
        legend('\mu = 5','\mu = 1','\mu = 0.5','\mu = 0');
        figure
    end
    system = tf(wn^2,[1,2*mu*wn,wn^2]);
    impulse(system)
    hold on
end
legend('\mu = -0.5');

%% 2.4, 2.5
wn = 10;
u = 0.5;
for Tp = logspace(-2,2,5)
    % implementation of the system
    TransferFcn = tf(1,[Tp,1]);
    TransferFcn1 = tf(wn^2,[1,2*u*wn,0]);
    Transfer = series(TransferFcn, TransferFcn1);
    unity = tf(1,1);
    system = feedback(Transfer, unity, -1);
    figure
    % bode plot
    subplot(2,2,[2,4])
    bode(system)
    title(['bode plot for T_p = ',num2str(Tp)]);
    % impulse response
    subplot(2,2,1)
    impulse(system)
    title(['impulse response for T_p = ',num2str(Tp)]);
    % pole-zero plot
    subplot(2,2,3)
    pzmap(system)
    title(['pole-zero plot for T_p = ',num2str(Tp)]);
end

%% 2.9
wn = 10;
u = 0.5;
for Tz = logspace(-2,2,5)
    % implementation of the system
    TransferFcn = tf([Tz,1],1);
    TransferFcn1 = tf(wn^2,[1,2*u*wn,wn^2]);
    system = series(TransferFcn, TransferFcn1);
    figure
    % bode plot
    subplot(2,2,[2,4])
    bode(system)
    title(['bode plot for T_z = ',num2str(Tz)]);
    % impulse response
    subplot(2,2,1)
    impulse(system)
    title(['impulse response for T_z = ',num2str(Tz)]);
    % pole-zero plot
    subplot(2,2,3)
    pzmap(system)
    title(['pole-zero plot for T_z = ',num2str(Tz)]);
end

%% Question03

Ts = 1;
table = [1, 2, 3];
for alpha = table
    system = tf([alpha], [1, alpha - 1], Ts);
    figure
    impulse(system)
    title(['$\alpha = T/RC = ', num2str(alpha), '$'], 'interpreter', 'latex');
    xlim([0, 20])
end

Ts = 1;
table = [0.1, 0.01];
for alpha = table
    system = tf([alpha], [1, alpha - 1], Ts);
    figure
    impulse(system)
    title(['$\alpha = T/RC = ', num2str(alpha), '$'], 'interpreter', 'latex');
    
end

%% Causal Two Point
Ts = 1;
for alpha = [0.01, 0.1, 1, 2, 3]
    system = tf([alpha, 0], [alpha + 1, -1], Ts);
    figure
    impulse(system)
    title(['$\alpha = T/RC = ', num2str(alpha), '$'], 'interpreter', 'latex');
end

%% Two point midpoint
alpha = 0.1;
Ts = 1;
system = tf([2*alpha, 0], [1,2*alpha, -1], Ts);
figure
impulse(system)
title(['$\alpha = T/RC = ', num2str(alpha), '$'], 'interpreter', 'latex');


%% Bilinear Mapping

for alpha = [0.01, 0.1, 1, 2, 3]
    Ts = 1;
    system = tf([alpha, alpha], [2+alpha, alpha-2], Ts);
    figure
    impulse(system)
    title(['$\alpha = T/RC = ', num2str(alpha), '$'], 'interpreter', 'latex');
end