clc
clear all
close all



%% verify c(z,t)
figure(101)
hold on
axis([0 1 0 1])
box on
xlabel('z')
ylabel('c, x')

figure(201)
hold on
axis([0 6 0 1])
box on
xlabel('z')
ylabel('c, x')
% plot([0 0.7], params.a0*[1 1], 'k--')

t = [linspace(1E-2, 0.5, 11), 2];
for i = 1:numel(t)
    params.t = t(i);
    params.a = sqrt(2*params.t);
    
    params.z2 = z2(params);
    
    if(t(i) <= 0.5)
        z = linspace(0, params.z2, 201);
        c = z/params.z2;
        x = (1-c)*params.z2;
        figure(101)
        plot(z,c, 'r-', 'LineWidth', 1)
        plot(z,x, 'k-', 'LineWidth', 1)
    else        
        z = params.z2 + linspace(0, 1, 3);        
        c = z - params.z2;
        x = 1-(z-params.z2);
        figure(201)
        plot([0,z,8],[0,c,1], 'r-', 'LineWidth', 1)
        plot([0,z,8],[1,x,0], 'k-', 'LineWidth', 1)
    end
end
return

%% verify G(x)/x
figure(300)
hold on
axis([0 2 0 params.g0])
box on
xlabel('x')
ylabel('G(x)/x')

params.t = 2;
params.a = sqrt(2*params.t);

x = x_grid(params, 0.0001);
out = G2X(x, params);

plot(x, out, 'k-', 'LineWidth', 1)
plot(params.a0*[1 1], [0 params.g0], 'k--')
plot([1 1], [0 params.g0], 'k--')

