clear all
close all 
clc
format long
%%
%parameters
alpha_deg = 8;%angle of attack [deg]
alpha = alpha_deg*2*pi/360;%angle of attack [rad]
theta = pi/6;%max delta wing inclination [rad]
U_inf = 60;%free stream velocity [m/s]
k_inf = pi;%shape coefficient [-]
rho = 1.2;%air density [kg/m^3]
L = 4;%mid span length [m]
lambda = 8;%aspect ratio [-]
S = 4*(L^2)/lambda;%projected surface for the rectangular chord distribution [m^2]
%%
%initialisation
N_max = 2e3;%maximum number of sampling points [-]

for wing = 1:3
    %%
    %chord length distribution along the span
    switch wing%if wing == 1 --> rectangular wing; if wing == 2 --> delta wing; if wing == 3 --> elliptic wing
        
        case 1
        
        c0 = 2*L/lambda;%chord length at y0=0 for the rectangular wing [m]
        l=@(y)(c0.*ones(size(y)));
        string_tit = 'Rectangular chord length distribution $l(y)$ $[m]$';
        
        case 2

            if atan(c0/L) > theta

                error('angle of delta wing higher than 30 degrees')

            end

            c0 = S/L;%chord length at y0=0 for the triangular wing [m]
            l = @(y)(-abs(c0.*y./L)+c0);
            string_tit = 'Triangular chord length distribution $l(y)$ $[m]$';

        case 3

            c0 = 2*S/(pi*L);%chord length at y0=0 for the triangular wing [m]
            l = @(y)(c0.*sqrt(1-(y./L).^2));
            string_tit = 'Elliptic chord length distribution $l(y)$ $[m]$';
            
    end
    
    err_rel_gamma = zeros(N_max-2,1);
    err_rel_lift = zeros(N_max-2,1);
    err_rel_drag = zeros(N_max-2,1);
    for N = 2: N_max
        y = linspace(-L,L,N)';%discretized y-direction
        eta = linspace(-L+(y(2)-y(1))*0.5,L-(y(end)-y(end-1))*0.5,N-1);%discrete control points
        [gamma w] = gamma_downwash_comp(y,eta,l,U_inf,alpha);%circulation along the span [m^2/s]
        gamma_tot = trapz(y,gamma);
        lift = rho*U_inf.*trapz(y,gamma);%lift force distribution along the span [N/m]
        drag = -2*trapz(y,w.*gamma)/(S*U_inf^2);%drag force distribution along the span [N/m]
        if N >= 3
            err_rel_gamma(N-2) = (gamma_tot - gamma_k)/gamma_k;
            err_rel_lift(N-2) = abs((lift - lift_k)/lift_k);
            err_rel_drag(N-2) = abs((drag - drag_k)/drag_k);
        end
        gamma_k = gamma_tot;
        lift_k = lift;
        drag_k = drag;
               
    end
    
    conv_values = ['\n \n',string_tit,'\n \n relative error for circulation: ',num2str(err_rel_gamma(end)),'\n \n relative error for lift: ',num2str(err_rel_lift(end)),'\n \n relative error for drag: ',num2str(err_rel_lift(end)),'\n'];
    fprintf(conv_values)
    %%
    %convergence plots for each chord length distribution
    figure
    grid on
    hold on
    tit = ['Total circulation $\Gamma$ $[\frac{m^{2}}{s}]$ qualitative convergence analysis for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    subtitle(string_tit,'interpreter','latex')
    xlabel('Number of sampling points along the span $N$ $[-]$','interpreter','latex')
    ylabel('$\frac{|\Gamma(y)_{k+1}-\Gamma(y)_{k}|}{|\Gamma(y)_{k}|}$ $[-]$','FontSize',18,'interpreter','latex')
    plot((3:1:N_max)',err_rel_gamma,'r^--')
    axis([0 100 0 max(err_rel_gamma)])
    hold off
    
    figure
    grid on
    hold on
    tit = ['Total lift force $L$ $[N]$ qualitative convergence analysis for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    subtitle(string_tit,'interpreter','latex')
    xlabel('Number of sampling points along the span $N$ $[-]$','interpreter','latex')
    ylabel('$\frac{|L_{k+1}-L_{k}|}{|L_{k}|}$ $[-]$','FontSize',18,'interpreter','latex')
    plot((3:1:N_max)',err_rel_lift,'b^--')
    axis([0 100 0 max(err_rel_lift)])
    hold off
    
    figure
    grid on
    hold on
    tit = ['Total drag force $D$ $[N]$ qualitative convergence analysis for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    subtitle(string_tit,'interpreter','latex')
    xlabel('Number of sampling points along the span $N$ $[-]$','interpreter','latex')
    ylabel('$\frac{|D_{k+1}-D_{k}|}{|D_{k}|}$ $[-]$','FontSize',18,'interpreter','latex')
    plot((3:1:N_max)',err_rel_lift,'k^--')
    axis([0 100 0 max(err_rel_drag)])
    hold off
    
end

angles_attack_deg = (0:1:15)';%angle of attack [deg]
Cl = zeros(size(angles_attack_deg));
Cd = zeros(size(angles_attack_deg));
perf = zeros(size(angles_attack_deg));
for p = 1:length(angles_attack_deg)
    
    alpha_deg = angles_attack_deg(p);
    alpha = alpha_deg*2*pi/360;%angle of attack [rad]
    [gamma w] = gamma_downwash_comp(y,eta,l,U_inf,alpha);%circulation along the span [m^2/s]
    %%
    %lift coefficient computation
    Cl(p) = 2*trapz(y,gamma)/(S*U_inf); %lift coefficient [-]
    %%
    %induced drag coefficient computation
    Cd(p) = -2*trapz(y,w.*gamma)/(S*U_inf^2); %induced drag coefficient [-]
    %%
    %performances computation
    if Cd(p) > 0.5*eps
        perf(p) = Cl(p)/Cd(p); %wing performance [-]
    end
end

%%
%performances plot for a given chord length distribution and freestream velocity

figure
hold on
grid on
xlabel('$C_{L}$ $[-]$','interpreter','latex')
ylabel('$\frac{C_{L}}{C_{D}}$ $[-]$','interpreter','latex')
tit = ['Wing performances comparison $\frac{C_{L}}{C_{D}}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
title(tit,'interpreter','latex')
subtitle(string_tit,'interpreter','latex')
plot(Cl,perf,'b^--')
plot(Cl,pi*lambda./Cl,'r--')
axis([Cl(2) max(Cl) 0 max(perf)])
legend('\Big($\frac{C_{L}}{C_{D}}\Big)_{numerical}$ $[-]$','\Big($\frac{C_{L}}{C_{D}}\Big)_{theoretical}$ $[-]$','interpreter','latex')
hold off

figure
hold on
grid on
xlabel('$\alpha$ $[^{\circ}]$','interpreter','latex')
ylabel('$C_{L}$ $[-]$','interpreter','latex')
tit = ['Lift coefficients comparison $C_{L}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
title(tit,'interpreter','latex')
subtitle(string_tit,'interpreter','latex')
plot(angles_attack_deg,Cl,'b^--')
plot(angles_attack_deg,(2*k_inf/(1 + (2*k_inf/(pi*lambda))))*(2*pi/360).*angles_attack_deg,'r--')
axis([0 max(angles_attack_deg) 0 max(Cl)])
legend('$(C_{L})_{numerical}$ $[-]$','$(C_{L})_{theoretical}$ $[-]$','interpreter','latex')
hold off