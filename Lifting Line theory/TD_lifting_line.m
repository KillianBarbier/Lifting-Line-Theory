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
lambda = 8;%aspect ratio [-]
S = 4*L^2/lambda;%projected surface for the rectangular chord distribution [m^2]
%%
%initialisation
N = 1e3;%number of sampling points [-]
y = linspace(-L,L,N)';%discretized y-direction
eta = linspace(-L+(y(2)-y(1))*0.5,L-(y(end)-y(end-1))*0.5,N-1);%discrete control points

for wing = 1:3%if wing == 1 --> rectangular wing; if wing == 2 --> delta wing; if wing == 3 --> elliptic wing
    
    %%
    %chord length distribution along the span

    switch wing

        case 1

            c0 = 0.5*S/L;%chord length at y0=0 for the rectangular wing [m]
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

    Cl = zeros(size(angles_attack_deg));
    Cd = zeros(size(angles_attack_deg));
    perf = zeros(size(angles_attack_deg));

    for p = 1:length(angles_attack_deg)

        alpha_deg = angles_attack_deg(p);
        alpha = alpha_deg*2*pi/360;%angle of attack [rad]
        [gamma w] = gamma_downwash_comp(y,eta,l,U_inf,alpha);%circulation along the span [m^2/s]
        lift = rho*U_inf.*gamma;%lift force distribution along the span [N/m]
        %%
        %plot of the results

        if alpha_deg == 8

            figure
            grid on
            hold on
            tit = ['Circulation distribution along the span $\Gamma(y_{0})$ $[\frac{m^{2}}{s}]$ for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
            title(tit,'interpreter','latex')
            subtitle(string_tit,'interpreter','latex')
            xlabel('$y$ $[m]$','interpreter','latex')
            ylabel('$\Gamma(y_{0})$ $[\frac{m^{2}}{s}]$, $l(y)$ $[m]$','interpreter','latex')
            plot(y,l(y),'r--')
            plot(y,gamma,'b-.','Markersize',2)
            legend('$l(y)$ $[m]$','$\Gamma(y)$ $[\frac{m^{2}}{s}]$','interpreter','latex')
            hold off

            figure
            grid on
            hold on
            tit = ['Lift force distribution along the span $L(y_{0})$ $[N]$ for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
            title(tit,'interpreter','latex')
            subtitle(string_tit,'interpreter','latex')
            xlabel('$y$ $[m]$','interpreter','latex')
            ylabel('$L(y_{0})$ $[N]$, $l(y)$ $[m]$','interpreter','latex')
            plot(y,l(y),'r--')
            plot(y,lift,'c-.')
            legend('$l(y)$ $[m]$','$L(y)$ $[\frac{N}{m}]$','interpreter','latex')
            hold off

        end
        %%
        %lift coefficient computation
        Cl(p) = 2*trapz(y,gamma)/(S*U_inf); %lift coefficient [-]

    %     drag = -2*w.*gamma./(S*U_inf^2);%drag force distribution along the span [N/m]

        if alpha_deg == 8

            figure
            hold on
            grid on
            tit = ['Induced velocity along the span $w(y_{0})$ $[\frac{m}{s}]$ for $\alpha$ $=$ ',num2str(alpha_deg),'$^{\circ}$ and $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];        
            title(tit,'interpreter','latex');
            subtitle(string_tit,'interpreter','latex')
            xlabel('$y$ $[m]$','interpreter','latex')
            ylabel('$w(y_{0})$ $[\frac{m}{s}]$, $l(y)$ $[m]$','interpreter','latex')
            plot(y,l(y),'r--')
            plot(y,w,'k-.','Markersize',2)
            legend('$l(y)$ $[m]$','$w(y)$ $[\frac{m}{s}]$','interpreter','latex')
            hold off

        end
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
    %polar plot for a given chord length distribution and freestream velocity

    figure
    hold on
    grid on
    xlabel('$C_{D}$ $[-]$','interpreter','latex')
    ylabel('$C_{L}$ $[-]$','interpreter','latex')
    tit = ['Polar plot $C_{L}(C_{D})$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    subtitle(string_tit,'interpreter','latex')
    % axis([0 max(Cd) 0 max(Cd)])
    plot(Cd,Cl,'g-.')
    hold off
    %%
    %performances plot for a given chord length distribution and freestream velocity

    figure
    hold on
    grid on
    xlabel('$C_{L}$ $[-]$','interpreter','latex')
    ylabel('$\frac{C_{L}}{C_{D}}$ $[-]$','FontSize',16,'interpreter','latex')
    tit = ['Wing performances $\frac{C_{L}}{C_{D}}$ for $U_{\infty}$ $=$ ',num2str(U_inf),' $[\frac{m}{s}]$'];
    title(tit,'interpreter','latex')
    subtitle(string_tit,'interpreter','latex')
    axis([Cl(2) max(Cl) 0 max(perf)])
    plot(Cl,perf,'c-.')
    hold off

end
