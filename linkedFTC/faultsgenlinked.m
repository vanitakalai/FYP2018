%generates all possible fault matrix delta for linked sensors
function D_linked = faultsgenlinked(m,l,l_array)
%m=total number of sensors
%l=number of blocks %l_array contains size of each block

D_linked=zeros(m,m,2^l);

for i=1:2^l
        tmp=de2bi(i-1,l,'left-msb')-eps;
        D_linked(:,:,i)=diag(repelem(tmp,l_array));
end

end