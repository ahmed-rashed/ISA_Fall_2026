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

h_G_col=nan(7*N_layer,1);
T_col=nan(7*N_layer,1);
p_col=nan(7*N_layer,1);

% Gradient layers
for n=1:2:7
    n_grad_layer=(n+1)./2;
    ind_row=(n-1)*N_layer+1:n*N_layer;
    h_G_col(ind_row)=linspace(h_G0_row(n),h_G0_row(n+1),N_layer).';
    T_col(ind_row)=T0_row(n)+a0_vec(n_grad_layer).*(h_G_col(ind_row)-h_G0_row(n));
    p_col(ind_row)=p0_row(n).*(T_col(ind_row)./T0_row(n)).^-(g_0./a0_vec(n_grad_layer)/R);
end

% Isothermal layers
for n=2:2:6
    ind_row=(n-1)*N_layer+1:n*N_layer;
    h_G_col(ind_row)=linspace(h_G0_row(n),h_G0_row(n+1),N_layer).';
    T_col(ind_row)=repmat(T0_row(n),N_layer,1);
    p_col(ind_row)=p0_row(n).*exp(-g_0.*(h_G_col(ind_row)-h_G0_row(n))./R./T0_row(n));
end
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