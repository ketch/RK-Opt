function x=initial_guess(s,p,class,starttype)
%function x=initial_guess(s,p,class,starttype)
%
% Set initial guess for RK coefficients
% Includes some good initial guesses for optimal SSP methods
if ~ischar(starttype)
    x=starttype;
else
x=[];
switch class
    case 'irk'
        % Implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s)); 
                x(s+1:2*s-1)=rand(1,s-1); 
                x(2*s)=1-sum(x(s+1:2*s-1));
                x(2*s+1:2*s+s^2)=rand(1,s^2);
                x(2*s+s^2+1)=-0.01;
            case 'smart'
                x(1:s)=(1:s)/s;
                x(s+1:2*s)=1/s;
                for i=1:s
                    x(2*s+s*(i-1)+1:2*s+s*(i-1)+i-1)=1/s;
                    x(2*s+s*(i-1)+i)=1/(2*s);
                end
                x(2*s+s^2+1)=-0.01;
                x=x.*(1+rand(size(x))/4);
        end

    case 'irk5'
        % Implicit SSP methods of order >=5 always have one row of A
        % equal to zero
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1));    %c's
                x(s:2*s-2)=rand(1,s-1);        %b's
                x(2*s-1)=1-sum(x(s:2*s-2)); 
                x(2*s:2*s-1+s*(s-1))=rand(1,s*(s-1));  %A's
                x(2*s+s*(s-1))=-0.01;            %r
            case 'smart'
                x(1:s-1)=(1:s-1)/s;
                x(s:2*s-1)=1/s;
                for i=2:s
                    x(2*s+s*(i-2):2*s+s*(i-2)+i-2)=1/s;
                    x(2*s-1+s*(i-2)+i-1)=1/(2*s);
                end
                x(2*s+s*(s-1))=-0.1;
                x=x.*(1+rand(size(x))/4);
        end

    case 'dirk'
        %Diagonally implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s));                     %c's
                x(s+1:2*s-1)=rand(1,s-1);                   %b's
                x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
                x(2*s+1:2*s+s*(s+1)/2)=rand(1,s*(s+1)/2);   %A's
                x(2*s+s*(s+1)/2+1)=-0.01;                   %r
            case 'smart'
                x(1:s)=(1:s)/s;                             %c's
                x(s+1:2*s)=1/s;                             %b's
                x(2*s+1:2*s+s*(s+1)/2)=1/s;                 %A's
                j=0;
                for i=1:s
                  j=j+i;
                  x(2*s+j)=1/(2*s);                          %Diagonal A's (A_ii)
                end
                x(2*s+s*(s+1)/2+1)=-(s-(p-2)+sqrt(s^2-(p-2)))/2;                   %r
                x=x.*(1+rand(size(x))/10);
        end

    case 'sspdirk5'
        % Implicit SSP methods of order >=5 always have one row of A
        % equal to zero
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1));    % c
                x(s:2*s-2)=rand(1,s-1)/s;      % b
                x(2*s-1)=1-sum(x(s:2*s-2));    %b(end)
                x(2*s:2*s-2+s*(s+1)/2)=rand(1,s*(s+1)/2-1)/3.;
                x(2*s+s*(s+1)/2-1)=-0.01;
            case 'smart'
                nbz=4;
                %modified 2nd order
                x(1:s-1)=(1:s-1)/s;
                x(s:2*s-1)=1/(s-nbz);  
                %Zero some b's:
                x(2*s-nbz:2*s-1)=0;
                x(s)=x(s)/2;
                x(2*s:2*s-2+s*(s+1)/2)=1/s;
                j=0;
                for i=2:s
                  j=j+i;
                  x(2*s+j-i)=1/(2*s);                          %First column A's (A_i1)
                  x(2*s+j-i+1)=3/(2*s);                        %Second column A's (A_i2)
                  x(2*s-1+j)=1/(2*s);                          %Diagonal A's (A_ii)
                end
                x(2*s+s*(s+1)/2-1)=-0.01;
                x=x.*(1+rand(size(x))/10);
        end

    case 'sdirk'
        % Singly diagonally implicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s)=sort(rand(1,s));                     %c's
                x(s+1:2*s-1)=rand(1,s-1)/s;                   %b's
                x(2*s)=1-sum(x(s+1:2*s-1));                 %last b
                x(2*s+1:2*s+1+s*(s-1)/2)=rand(1,s*(s-1)/2+1)/s;   %A's
                x(2*s+1)=rand/3/s;
                x(2*s+1+s*(s-1)/2+1)=-0.01;                   %r
        end

    case 'erk'
        %Explicit Runge-Kutta
        switch starttype
            case 'random'
                x(1:s-1)=sort(rand(1,s-1)-1/2); 
                x(s:2*s-2)=rand(1,s-1)-1/2; 
                x(2*s-1)=1-sum(x(s:2*s-2));
                x(2*s:2*s-1+s*(s-1)/2)=rand(1,s*(s-1)/2)-1/2;
                x(2*s+s*(s-1)/2)=-0.01;
            case 'smart'
                r=s-(p-3)-sqrt(s-(p-3));
                x(1:s-1)=(1:s-1)/r;
                x(s:2*s-2)=1/s;
                x(2*s-1)=1-sum(x(s:2*s-2));
                x(2*s:2*s-1+s*(s-1)/2)=1/r;
                x(2*s+s*(s-1)/2)=-r;
        end

    otherwise
        % Low-storage methods
        n=set_n(s,class);
        x=rand(1,n);
end

end
