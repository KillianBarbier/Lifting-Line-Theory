clear all
close all 
clc
format long
%%
%parameters
angles_attack_deg = (0:1:15)';%angle of attack [deg]
theta = pi/6;%max delta wing inclination [rad]
U_inf = 6;%free stream velocity [m/s]
rho = 1.2;%air density [kg/m^3]
L = 4;%mid span length [m]
asp_ratios = (8:4:36)';%aspect ratios [-]
%%
%initialisation
N = 5e2;%number of sampling points [-]
y = linspace(-L,L,N)';%discretized y-direction
eta = linspace(-L+(y(2)-y(1))*0.5,L-(y(end)-y(end-1))*0.5,N-1);%discrete control points

Cl_wing = zeros(length(angles_attack_deg),3);
Cd_wing = zeros(length(angles_attack_deg),3);

perf_lambda = zeros(length(asp_ratios),3);
Cd_lambda = zeros(length(asp_ratios),3);
Cl_lambda = zeros(length(asp_ratios),3);

for wing = 1:3%if wing == 1 --> rectangular wing; if wing == 2 --> delta wing; if wing == 3 --> elliptic wing
    
    leg = cell(size(asp_ratios));
    
    figure(wing)
    grid on
    hold on
    xlabel('$C_{L}$ $[-]$','interpreter','latex')
    ylabel('$\frac{C_{L}}{C_{D}}$ $[-]$','FontSize',16,'interpreter','latex')
    tit = ['Wing performances $\frac{C_{L}}{C_{D}}$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    
    for k = 1:length(asp_ratios)
        
        Cl = zeros(size(angles_attack_deg));
        Cd = zeros(size(angles_attack_deg));
        perf = zeros(size(angles_attack_deg));
        
        lambda = asp_ratios(k);
        
        S = 4*L^2/lambda;%projected surface for the rectangular chord distribution [m^2]
        
        %%
        %chord length distribution along the span
    
        switch wing

            case 1

                c0 = S/(2*L);%chord length at y0=0 for the rectangular wing [m]
                l=@(y)(c0.*ones(size(y)));
                string_tit = 'Rectangular chord length distribution $l(y)$ $[m]$';

            case 2

                c0 = S/L;%chord length at y0=0 for the triangular wing [m]

                if atan(c0/L)>theta

                    error('angle of delta wing higher than 30 degrees')

                end

                l = @(y)(-abs(c0.*y./L)+c0);
                string_tit = 'Triangular chord length distribution $l(y)$ $[m]$';

            case 3

                c0 = 2*S/(pi*L);%chord length at y0=0 for the triangular wing [m]
                l = @(y)(c0.*sqrt(1-(y./L).^2));
                string_tit = 'Elliptic chord length distribution $l(y)$ $[m]$';

        end

        for p = 1:length(angles_attack_deg)

            alpha_deg = angles_attack_deg(p);
            alpha = alpha_deg*2*pi/360;%angle of attack [rad]
            [gamma w] = gamma_downwash_comp(y,eta,l,U_inf,alpha);%circulation along the span [m^2/s]
            lift = rho*U_inf.*gamma;%lift force distribution along the span [N/m]
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
            
            if alpha_deg == 10;
                
                perf_lambda(k,wing) = perf(p);
                Cd_lambda(k,wing) = Cd(p);
                Cl_lambda(k,wing) = Cl(p);
                
                alpha_comp = alpha_deg;
                
            end

        end

    leg{k} = ['$\lambda$ $=$',num2str(asp_ratios(k)),' $[-]$'];
    plot(Cl,perf,'s-.')
    
    if lambda == asp_ratios(end)
        
        Cl_wing(:,wing) = Cl;
        Cd_wing(:,wing) = Cd;
        
    end
    
    end
    
    subtitle(string_tit,'interpreter','latex')
    axis([Cl(2) max(Cl) 0 max(perf)])
    legend(leg,'interpreter','latex')
    hold off
        
end

%%
%comparison of Cl vs alpha for the three types of wings

figure
hold on
grid on
xlabel('$\alpha$ $[^{\circ}]$','interpreter','latex')
ylabel('$C_{L}$ $[-]$','interpreter','latex')
tit = ['Lift coefficient comparison $C_{L}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$ and $\lambda$ $=$ ',num2str(lambda)];
title(tit,'interpreter','latex')
plot(angles_attack_deg,Cl_wing(:,1),'--o')
plot(angles_attack_deg,Cl_wing(:,2),'--*')
plot(angles_attack_deg,Cl_wing(:,3),'--^')
legend('Rectangular chord distribution $l(y)$','Triangular chord distribution $l(y)$','Elliptic chord distribution $l(y)$','Location','northwest','interpreter','latex')
hold off

%%
%comparison of Cd vs alpha for the three types of wings

figure
hold on
grid on
xlabel('$\alpha$ $[^{\circ}]$','interpreter','latex')
ylabel('$C_{D}$ $[-]$','interpreter','latex')
tit = ['Drag coefficient comparison $C_{D}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$ and $\lambda$ $=$ ',num2str(lambda)];
title(tit,'interpreter','latex')
plot(angles_attack_deg,Cd_wing(:,1),'--o')
plot(angles_attack_deg,Cd_wing(:,2),'--*')
plot(angles_attack_deg,Cd_wing(:,3),'--^')
legend('Rectangular chord distribution $l(y)$','Triangular chord distribution $l(y)$','Elliptic chord distribution $l(y)$','Location','northwest','interpreter','latex')
hold off

%%
%Comparison of performances for the three type of wing as a function of
%lambda

figure
hold on
grid on
xlabel('$\lambda$ $[-]$','interpreter','latex')
ylabel('$\frac{C_{L}}{C_{D}}$ $[-]$','FontSize',16,'interpreter','latex')
tit = ['Performances comparison $\frac{C_{L}}{C_{D}}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$ and $\alpha$ $=$ ',num2str(alpha_comp),' $^{\circ}$'];
title(tit,'interpreter','latex')
plot(asp_ratios,perf_lambda(:,1),'--o')
plot(asp_ratios,perf_lambda(:,2),'--*')
plot(asp_ratios,perf_lambda(:,3),'--^')
axis([min(asp_ratios) max(asp_ratios) 0 max([perf_lambda(:,1);perf_lambda(:,2);perf_lambda(:,3)])])
legend('Rectangular chord distribution $l(y)$','Triangular chord distribution $l(y)$','Elliptic chord distribution $l(y)$','Location','northwest','interpreter','latex')
hold off

%%
%Comparison of Cd for the three type of wing as a function of
%lambda

figure
hold on
grid on
xlabel('$\lambda$ $[-]$','interpreter','latex')
ylabel('$C_{D}$ $[-]$','interpreter','latex')
tit = ['Drag coefficient comparison $C_{D}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$ and $\alpha$ $=$ ',num2str(alpha_comp),' $^{\circ}$'];
title(tit,'interpreter','latex')
plot(asp_ratios,Cd_lambda(:,1),'--o')
plot(asp_ratios,Cd_lambda(:,2),'--*')
plot(asp_ratios,Cd_lambda(:,3),'--^')
axis([min(asp_ratios) max(asp_ratios) 0 max([Cd_lambda(:,1);Cd_lambda(:,2);Cd_lambda(:,3)])])
legend('Rectangular chord distribution $l(y)$','Triangular chord distribution $l(y)$','Elliptic chord distribution $l(y)$','Location','northeast','interpreter','latex')
hold off

%%
%Comparison of Cd for the three type of wing as a function of
%lambda

figure
hold on
grid on
xlabel('$\lambda$ $[-]$','interpreter','latex')
ylabel('$C_{L}$ $[-]$','interpreter','latex')
tit = ['Lift coefficient comparison $C_{D}$ $[-]$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$ and $\alpha$ $=$ ',num2str(alpha_comp),' $^{\circ}$'];
title(tit,'interpreter','latex')
plot(asp_ratios,Cl_lambda(:,1),'--o')
plot(asp_ratios,Cl_lambda(:,2),'--*')
plot(asp_ratios,Cl_lambda(:,3),'--^')
axis([min(asp_ratios) max(asp_ratios) 0 max([Cl_lambda(:,1);Cl_lambda(:,2);Cl_lambda(:,3)])])
legend('Rectangular chord distribution $l(y)$','Triangular chord distribution $l(y)$','Elliptic chord distribution $l(y)$','Location','southeast','interpreter','latex')
hold off