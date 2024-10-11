function label = make_label(outG0)
[n, K] = size(outG0); label = zeros(n,1); 
for i = 1:n
    for j = 1:K
        if outG0(i,j)==1
            label(i) = j;
        end
    end
end

