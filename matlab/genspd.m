function [S, invS] = genspd(d, density, m, v)

% GENSPD Generate a random sparse positive-definite matrix.
%
% FORMAT
% DESC Generates a random sparse positive-definite matrix w/ given density and
% non-zero entries normally distributed w/ given mean and variance.
% ARG d : size of the sparse matrix.
% ARG density : density of the sparse matrix; the ratio of non-zero entries to all
% entries, above the diagonal.
% RETURN S : the generated sparse positive-definite matrix.
% RETURN Sinv : the inverse of the generate 
%
% SEEALSO : chol, pdinv
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

spN = ceil((d^2-d)/2 * density);    % Number of non-zeros elements that satisfy the density.
sel = randperm(d^2);                % Random positions.
sel( mod(sel,d+1) == 1 ) = [];      % Remove positions of the diagonal.
sel = sel( 1:spN );                 % Enough random positions to satisfy the density.
S = zeros(d);
S(sel) = 1;                         % Non-zero positions mask.
S = triu(S + S',1);                 % Upper triangular after balancing.
S = S .* (randn(d)*sqrt(v) + m);    % Apply mask to random non-zero entries.
S = S + S';                         % Symmetrify.
pd = 1;
diagonal = eye(d);
while pd > 0
    testS = S + diagonal;
    [C, pd] = chol(testS);
    diagonal = diagonal*2;
end
S = testS;

if nargout > 1
    invC = C\eye(d);
    invS = invC*invC';
end


