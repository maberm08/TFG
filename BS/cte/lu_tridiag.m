function [subdiag_L, diag_U, superdiag_U] = lu_tridiag(subdiag_A, diag_A, superdiag_A)
% LU_TRIDIAG
% Factorización LU de una matriz tridiagonal sin pivoteo (O(n)).
%
% ENTRADAS:
%   subdiag_A   : subdiagonal de A, de longitud n-1
%   diag_A      : diagonal principal de A, de longitud n
%   superdiag_A : superdiagonal de A, de longitud n-1
%
% SALIDAS:
%   subdiag_L   : subdiagonal de L, con diag(L)=1
%   diag_U      : diagonal principal de U
%   superdiag_U : superdiagonal de U

    n = length(diag_A);

    subdiag_L   = zeros(n-1,1);
    diag_U      = zeros(n,1);
    superdiag_U = superdiag_A;

    diag_U(1) = diag_A(1);

    for k = 1:n-1
        subdiag_L(k) = subdiag_A(k) / diag_U(k);
        diag_U(k+1)  = diag_A(k+1) - subdiag_L(k) * superdiag_U(k);
    end
end