function [A, expanded_vector, expanded_matrix] = ...
    normalise_matrix(A,direction, expanded_vector, expanded_matrix)
%This function normalises along each row of a matrix

if(strcmpi(direction,'r'))
    %expanded_matrix = abs(A).*abs(A);
    expanded_matrix = bsxfun(@times,abs(A),abs(A));
    expanded_vector = sqrt(sum(expanded_matrix,2)); %norm is square root of sum of squares for each row
    A = bsxfun(@rdivide,A,expanded_vector); %Divide each value in row by norm of value in that row

elseif(strcmpi(direction,'c'))
    %expanded_matrix = abs(A).*abs(A);
    expanded_matrix = bsxfun(@times,abs(A),abs(A));
    expanded_vector = sqrt(sum(expanded_matrix,1)); %norm is square root of sum of squares for each column
    A = bsxfun(@rdivide,A,expanded_vector); %Divide each value in column by norm of value in that column

else
    error('Choose ''r'' for normalising across rows or ''c'' for down columns');
end


end