ncols=length(cols);
B = sym(zeros(p+1,cols));
d = sym(zeros(p+1,1));

M=k*(s+1);
count=1;
for i=0:p
    d(i+1)=k^i;
    for m=1:M
        if m==cols(count)
            ii=ceil(m/(s+1));
            j=m-(ii-1)*(s+1)-1;
            B(i+1,count)=0;
            for l=0:min(i,j)
                B(i+1,count)=B(i+1,count)+factorial(i)/factorial(l)/factorial(i-l) ...
                  *(k-ii)^(i-l)/r^l*prod(j-(0:(l-1)));
            end
        count = count+1;
        end
    end
end

X = A\b
