function Y = getY(y,m)

idx = 0;
Y = zeros(m);
for col = 1 : m
    for row = col : m
        idx = idx + 1;
        Y(row, col) = y(idx);
    end
end

Y = Y + Y' - diag(diag(Y));