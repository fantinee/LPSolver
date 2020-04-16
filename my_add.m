function [c] = my_add(a,b)
c = a+my_mult(a,b);
end

function [c] = my_mult(a,b)
c = a*b;
end
