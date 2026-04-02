
fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)

clc
close all
clear all

%%
m = 100701;
eps = 1e-7;
for n = [0,4,9, 15]
    y = linspace(-1+eps, 1-eps, m);
    eigvals = zeros(2, m);
    for i = 1:m
        eigvals(:, i) = ...
            (eig(Jac_singular_legendre(y(i), [], n)));
    end
    
    figure(n+1)
    hold on
    box on
    grid on
    axis([-1 1 -150 150])
    h1 = ...
        plot(y, imag(eigvals(1,:)), 'r.', 'DisplayName', 'imag');
    h2 = ...
        plot(y, real(eigvals(1,:)), 'k.', 'DisplayName', 'real');
    
    plot(y, imag(eigvals), 'r.')
    plot(y, real(eigvals), 'k.')
    legend([h1, h2])
    xlabel('{\itx}')
    ylabel('eigenvals')
    title(['{\lambda} = {\itn}({\itn}+1), {\itn} = ', num2str(n)])
    saveas(n+1, ['Figs\eig_', num2str(n), '.emf'])
end
