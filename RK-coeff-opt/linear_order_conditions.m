function linear_cond = linear_order_conditions(A,b,p)
%linear_cond = linear_order_conditions(A,b,p)
%
% Check linear order conditions

linear_cond = b'*A^(p-1)*ones(length(b),1) - 1/factorial(p);

end

