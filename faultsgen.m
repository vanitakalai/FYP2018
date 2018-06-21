%generates all possible fault matrix delta
function D = faultsgen(m)

D=zeros(m,m,2^m);
for i=1:2^m
    D(:,:,i)=diag(de2bi(i-1,m,'left-msb')-eps);
end

end
