function save_figs(sigma, params)
fig_id = 700;

for i = 1:6
    figure(fig_id+i)
    title(['{\sigma}=', num2str(sigma), ', f=', num2str(params.f), ', R=' num2str(params.R)])
end

exp = {'emf', 'fig'};
vars = {'Phi', 'Psi', 'Chi', 'Gam', 'Omega', 'Y'};

for i = 1:2
    e = exp{i};
    for j = 1:6
        v = vars{j};
        saveas(fig_id+j, ...
            ['Figs/',e,'/',v,'_sigma=', caption(sigma, params), ['.' e]]);
    end
end

end


function out = caption(sigma, params)
out = [num2str(sigma), '_f=', num2str(params.f), '_R=' num2str(params.R)];
end