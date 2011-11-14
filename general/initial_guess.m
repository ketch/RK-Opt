function x=initial_guess_ussp(s,class,starttype)
%function x=initial_guess_ussp(s,class,starttype)
%
% Set initial guess for RK coefficients
%This function needs to be cleaned up!
x=[];
switch class

  case 'irk'
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
        x(2*s+s^2+1)=-0.01;%(s-(p-2)+sqrt(s^2-(p-2)));
	x=x.*(1+rand(size(x))/4);
    end

  case 'irk5'
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
      case 'load'
	load dirk66b A b c r;
	x=loadX(A,b,c,r,'irk5');
    end

  case 'dirk'
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
%	load dirk114.mat A b c r;
%	x=loadX(A,b,c,r,'dirk');
%	x=x.*(1+(rand(size(x))-1/2)/8);
      case 'load'
	x=loadX(A,b,c,r,'dirk');
    end

  case 'dirk5'
    switch starttype
      case 'mod'  %Load another method and modify it
        load optimal_RK_schemes/implicit/dirk45_newbest.mat A b c r;
	A=[A;rand(1,s-1)/50.]; A=[A zeros(s,1)];
	A(end,end)=rand/50.;
	b=[b rand/50]; c=[c rand/50];
	x=Abc_to_x(A,b,c,r,'dirk5');
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

	x(s)=x(s)/2;% x(s+1)=x(s+1)*2;
        x(2*s:2*s-2+s*(s+1)/2)=1/s;
	j=0;
	for i=2:s
          j=j+i;
%          x(2*s+j-i+i-nbz:2*s+j-i+i-1)=0;
          x(2*s+j-i)=1/(2*s);                          %First column A's (A_i1)
          x(2*s+j-i+1)=3/(2*s);                        %Second column A's (A_i2)
          x(2*s-1+j)=1/(2*s);                          %Diagonal A's (A_ii)
	end
        %modified 3rd order
%        x(1:s-1)=0.5*(1-sqrt((s-1)/(s+1)))+(1:s-1)/sqrt(s^2-1);
%	x(s:2*s-1)=1/s;
%        x(2*s:2*s-2+s*(s+1)/2)=1/sqrt(s^2-1);
%	j=0;
%	for i=2:s
%          j=j+i;
%          x(2*s-1+j)=0.5*(1-sqrt((s-1)/(s+1)));
%          x(2*s+j-i)=1/(3*s);                          %First column A's (A_i1)
%	end
        x(2*s+s*(s+1)/2-1)=-0.01;
	x=x.*(1+rand(size(x))/10);
      case 'load'
	load bestdirk76.mat A b c r;
	x=Abc_to_x(A,b,c,r,'dirk5');
%	x=x.*(1+(rand(size(x))-1/2)/100);
    end

  case 'sdirk'
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
    switch starttype
      case 'random'
        x(1:s-1)=sort(rand(1,s-1)-1/2); 
	x(s:2*s-2)=rand(1,s-1)-1/2; 
	x(2*s-1)=1-sum(x(s:2*s-2));
        x(2*s:2*s-1+s*(s-1)/2)=rand(1,s*(s-1)/2)-1/2;
%	x=-x;
        x(2*s+s*(s-1)/2)=-0.01;
      case 'smart'
        r=s-(p-3)-sqrt(s-(p-3));
        x(1:s-1)=(1:s-1)/r;
	x(s:2*s-2)=1/s;
	x(2*s-1)=1-sum(x(s:2*s-2));
        x(2*s:2*s-1+s*(s-1)/2)=1/r;
        x(2*s+s*(s-1)/2)=-r;
%	load optimal_RK_schemes/erk/erk174.mat A b c r;
	load erk184col6.mat A b c r;
	x=loadX(A,b,c,r,'erk');
%	x=x.*(1+(rand(size(x))-1/2)/100);
%	x=x.*(1+(rand(size(x))-1/2)/5);
      end
    otherwise
      n=set_n(s,class);
      x=rand(1,n);
end


