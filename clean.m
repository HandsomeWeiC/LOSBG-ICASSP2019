function clean(Z,k)
[m,n]=size(Z);
S=Z;
if k==1
    for j=1:n
        for i=1:m
            if S(i,j)~=0
                S(i,j)=1;
            end
        end
    end
else
    for j=1:n
        for i=1:m
            a=10*S(i,j);
            if a<1
                S(i,j)=a;
            else
                S(i,j)=1;
            end
        end
    end
end

     
figure
% imshow(S)
imshow(S);
