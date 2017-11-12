function ind = maxind(x)
% % find local maxima in x
% % x should be a vector

if size(x,2) == 1
    x = x';
end

m = length(x);
x_rshift = [x(1) x(1:m-1)];
x_lshift = [x(2:m) x(m)];
x_max = (x > x_rshift) & (x > x_lshift);
ind = find(x_max == 1);

ind( find(ind==1) ) = [];
ind( find(ind==m) ) = [];