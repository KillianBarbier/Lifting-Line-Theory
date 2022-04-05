function [gamma w] = gamma_downwash_comp(y,eta,l,U_inf,alpha)

%%
%numerical resolution for the circulation along the span

%initialisation
N = length(y);
A = zeros(N,N);
b = zeros(N,1);

%Boundary Conditions setting
A(1,1) = 1;
A(N,N) = 1;
b(1) = 0;
b(N) = 0;

for i = 2:N-1
    b(i) = pi*U_inf*alpha*l(y(i));
    for j = 2:N-1
        if i == j
            A(i,j) = l(y(i))*0.25*(1/(eta(j)-y(i)) - 1/(eta(j-1)-y(i))) + 1;
        else
            A(i,j) = l(y(i))*0.25*(1/(eta(j)-y(i)) - 1/(eta(j-1)-y(i)));
        end
    end
end

gamma = A\b;%circulation along the span [m^2/s]

%%
%induced velocity computation

%initialisation
w = zeros(size(gamma));

for i = 1:N
    for j = 1:N-1
        w(i) = w(i) + 0.25*((gamma(j+1)-gamma(j))/(eta(j)-y(i)))/pi;
    end
end

end