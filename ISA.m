clearvars
clc
close all

g_0=9.80665;
R=287.04;
r=6.356766e6;
gamma=1.4;

h_G0_row=[0,11,25,47,53,79,90,105]*1e3;
T0_row=[288.16,216.66,216.66,282.66,282.66,165.66,165.66];
p0_row=[101330,22633,2488.7,120.45,58.323,1.0095,.10444];
a0_vec=[-.0065,.003,-.0045,.004];

N_layer=50;

% Gradient layers
h_G1_col=linspace(h_G0_row(1),h_G0_row(2),N_layer).';
T1_col=T0_row(1)+a0_vec(1).*(h_G1_col-h_G0_row(1));
p1_col=p0_row(1).*(T1_col./T0_row(1)).^-(g_0./a0_vec(1)/R);

h_G3_col=linspace(h_G0_row(3),h_G0_row(4),N_layer).';
T3_col=T0_row(3)+a0_vec(2).*(h_G3_col-h_G0_row(3));
p3_col=p0_row(3).*(T3_col./T0_row(3)).^-(g_0./a0_vec(2)/R);

h_G5_col=linspace(h_G0_row(5),h_G0_row(6),N_layer).';
T5_col=T0_row(5)+a0_vec(3).*(h_G5_col-h_G0_row(5));
p5_col=p0_row(5).*(T5_col./T0_row(5)).^-(g_0./a0_vec(3)/R);

h_G7_col=linspace(h_G0_row(7),h_G0_row(8),N_layer).';
T7_col=T0_row(7)+a0_vec(4).*(h_G7_col-h_G0_row(7));
p7_col=p0_row(7).*(T7_col./T0_row(7)).^-(g_0./a0_vec(4)/R);

% Isothermal layers
h_G2_col=linspace(h_G0_row(2),h_G0_row(3),N_layer).';
T2_col=repmat(T0_row(2),N_layer,1);
p2_col=p0_row(2).*exp(-g_0.*(h_G2_col-h_G0_row(2))./R./T0_row(2));

h_G4_col=linspace(h_G0_row(4),h_G0_row(5),N_layer).';
T4_col=repmat(T0_row(4),N_layer,1);
p4_col=p0_row(4).*exp(-g_0.*(h_G4_col-h_G0_row(4))./R./T0_row(4));

h_G6_col=linspace(h_G0_row(6),h_G0_row(7),N_layer).';
T6_col=repmat(T0_row(6),N_layer,1);
p6_col=p0_row(6).*exp(-g_0.*(h_G6_col-h_G0_row(6))./R./T0_row(6));

h_G_col=[h_G1_col;h_G2_col;h_G3_col;h_G4_col;h_G5_col;h_G6_col;h_G7_col];
T_col=[T1_col;T2_col;T3_col;T4_col;T5_col;T6_col;T7_col];
p_col=[p1_col;p2_col;p3_col;p4_col;p5_col;p6_col;p7_col];

h_col=r.*h_G_col./(r+h_G_col);
rho_col=p_col./R./T_col;
a_vec=sqrt(gamma.*R.*T_col);

figure;
plot(h_col/1e3,h_G_col/1e3)
xlabel('$h$ (km)','interpreter','latex') 
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(T_col-273,h_G_col/1e3)
xlim("padded")
grid on 
xlabel('$T$ (\r{ }C)','interpreter','latex') 
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(a_vec/1e3*60*60,h_G_col/1e3)
xlim("padded")
grid on 
xlabel('$a$ (km/h)','interpreter','latex') 
ylabel('$h_{G}$ (km)','interpreter','latex')

figure;
tiledlayout(1,2)
nexttile
plot(p_col/1e5,h_G_col/1e3)
xlim("padded")
grid on 
xlabel('$p$ (bar)','interpreter','latex') 
ylabel('$h_{G}$ (km)','interpreter','latex')

nexttile
plot(rho_col,h_G_col/1e3)
xlim("padded")
grid on 
xlabel('$\rho$ (kg/m\textsuperscript{3})','interpreter','latex') 
ylabel('$h_{G}$ (km)','interpreter','latex')