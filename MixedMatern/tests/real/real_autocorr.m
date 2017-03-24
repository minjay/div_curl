% check independence

load('wind.mat')
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.075 0.02], [0.075 0.02]);
u = samples(:, 1:2:end);
v = samples(:, 2:2:end);
rng(1)
index = randsample(1:n, 16);
for i = 1:16
    h = subplot(4, 4, i);
    autocorr(u(:, index(i)))
    delete(findall(h,'Type','text'))
    if i>=13
        xlabel('Lag')
    end
    if mod(i, 4)==1
        ylabel('Sample Autocorr')
    end
end

figure
for i = 1:16
    h = subplot(4, 4, i);
    autocorr(v(:, index(i)))
    delete(findall(h,'Type','text'))
    if i>=13
        xlabel('Lag')
    end
    if mod(i, 4)==1
        ylabel('Sample Autocorr')
    end
end


