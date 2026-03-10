function [A,b,c,A_2N,B_2N]=unpack_2N(X)
    % n = 2s - 1 free parameters
    s = (length(X)+1)/2; % # of stages
    % for 2N methods:
    % X = [A_2, A_3, ..., A_s, B_1, ..., B_s]
    A_2N = [0 X(1:s-1)];
    B_2N = X(s:end);
    K = zeros(s+1,s);
    for i=1:s-1
        K(i+1,i) = B_2N(i);
        for j = 2:s-i+1
            for k=0:j-1
                K(i+j,i) = K(i+j,i) + B_2N(i+k)*prod(A_2N(i+1:i+k));
            end
        end
    end
    A = K(1:s,:);
    b = K(s+1,:); b=b'; b(end)=B_2N(end);
    c=sum(A,2);
    alpha=zeros(s+1,s); beta=zeros(s+1,s); % Not yet implemented


