function qqchi2(y)
% Provide Chi-square QQ plot for multivariate data
%
% y: data matrix
%
[nr,nc]=size(y);
if nr <= 1
    stop
end
avey=mean(y);
dev = y - kron(ones(nr,1),avey);
s = dev'*dev/(nr-1);
si = inv(s);
d2=diag(dev*si*dev');
d3=sort(d2);
p=[];
for i=1:nr
    p=[p;i];
end
prob=(p-0.5)/nr;
q1 = chi2inv(prob,nc);
plot(q1,d3,'+');
hold on
plot([min(q1), max(q1)], [min(q1), max(q1)], 'r-.')
box off
