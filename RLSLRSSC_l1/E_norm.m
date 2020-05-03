function [value] = E_norm(E, error_type)
if error_type == 1
    value= norm(E,1);
else
    %l21_norm
    [d,N] = size(E);
    temp = 0;
    for i = 1:N
        i_norm = norm(E(:,i));
        temp = max(temp,i_norm);
    end
    value = temp;
end
end