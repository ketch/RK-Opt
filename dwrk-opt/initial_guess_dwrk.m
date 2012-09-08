function x=initial_guess_dwrk(s,class,starttype)
% function x=initial_guess_dwrk(s)
%
% Set initial guess
% Options:
% smart  - based on known good/optimal methods
% random - random among some promising region
% load   - restart from a method found before
% ::
%       x = [c' b' A ct' bt' At r] 
%
% (A,At are stored row-by-row)

[n,n2]=set_n_dwrk(s,class);
x=zeros(n,1);

switch class
  case 'erk'
    x(1:s-1)=sort(rand(1,s-1)-1/2);                % c
    x(1+n2:s-1+n2)=sort(rand(1,s-1)-1/2);          % ct
    x(s:2*s-2)=rand(1,s-1)-1/2;                    % b
    x(s+n2:2*s-2+n2)=rand(1,s-1)-1/2;              % bt
    x(2*s-1)=1-sum(x(s:2*s-2));                    % b_s
    x(2*s-1+n2)=-1-sum(x(s+n2:2*s-2+n2));           % bt_s
    x(2*s:2*s-1+s*(s-1)/2)=rand(1,s*(s-1)/2)-1/2;  % A
    x(2*s+n2:2*s-1+s*(s-1)/2+n2)=rand(1,s*(s-1)/2)-1/2;  %At
  case 'dirk'
    x = rand(n,1) - 1/2;
    x(1   :s   ) = sort(rand(1,s)-1/2);          % c
    x(1+n2:s+n2) = sort(rand(1,s)-1/2);          % ct
    x(s+1   :2*s   ) = rand(1,s)-1/2;            % b
    x(s+1+n2:2*s+n2) = rand(1,s)-1/2;            % bt
  case 'irk'
    if strcmp(starttype,'smart')
        assert(s==2,'Use s=2!');
        r=10.;
        x(1)   =(r^2-2*r-2.)/(2*r);                    % c_1
        x(2)   =(r-2.)/2.;                            % c_2
        x(1+n2)=(r^2-4*r+2.)/(2*r);                      % ct_1
        x(2+n2)=(r^2-4*r+2.)/(2*r);                      % ct_2
        x(s+1)=(r-2.)/2.;                          % b_1
        x(2*s)=1./r;                               % b_2
        x(s+n2+1)=0.;                              % bt_1
        x(s+n2+2)=(r^2-4*r+2.)/(2*r);              % bt_2
        x(2*s+1)=(r^2-2*r-2.)/(2*r);               % A_11
        x(2*s+2)=0;                                % A_12
        x(2*s+3)=(r-2.)/2;                         % A_21
        x(2*s+4)=0;                                % A_22
        x(n2+2*s+1)=0.;                            % At_11
        x(n2+2*s+1)=(r^2-4*r+2.)/(2*r);            % At_12
        x(n2+2*s+1)=0.;                            % At_21
        x(n2+2*s+1)=(r^2-4*r+2.)/(2*r);            % At_22
    elseif strcmp(starttype,'random')
        x(1:s)=sort(rand(1,s)-1/2);                % c
        x(1+n2:s+n2)=sort(rand(1,s)-1/2);          % ct
        x(s+1:2*s)=rand(1,s)-1/2;                  % b
        x(s+n2+1:2*s+n2)=rand(1,s)-1/2;            % bt
        x(2*s)=1-sum(x(s+1:2*s-1));                % b_s
        x(2*s+n2)=-1-sum(x(s+1+n2:2*s-1+n2));      % bt_s
        x(2*s+1:2*s+s^2)=rand(1,s^2)-1/2;           % A
        x(n2+2*s+1:n2+2*s+s^2)=rand(1,s^2)-1/2;     % At
    end

end

x(n)=-0.01;                                    % r

