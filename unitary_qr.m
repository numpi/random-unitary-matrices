function [d, its] = unitary_qr(G, d)
%UNITARY_QR Compute the eigenvalues of a unitary matrix. 
%
% D = UNITARY_QR(G, D) computes the eigenvalues of the unitary matrix
%
%   Q = G_1 ... G_{n-1} diag(D), 
%
% where G_j is a matrix differing from the identity only on the diagonal
% block of indice (j, j+1), which is equal to G{j}, and diag(D) is diagonal 
% with entries of modulus one given by the vector D. 
%
% The core factors G{j} need to be unitary themselves. 

if exist('unitary_qr_eiscor', 'file') ~= 3
    fprintf('UNITARY_QR: Cannot find the MEX interface to eiscor\n |\n');
    fprintf(' | To compile it, you can download the EISCOR library and \n');
    fprintf(' | compile it by running the compile_eiscor script.\n\n');
    
    return;
end

% Given a sequence of Givens rotations G{1} ... G{n-1}, we compute the
% eigenvalues of their product.
n = length(G) + 1;

if ~exist('d', 'var')
    d = ones(n, 1);
end

% Make sure we have rotations with real sines
[G, d] = normalize_rotations(G, d);

% Call eiscor, for which we need to cast everything back to real data. 
Q = zeros(3, n - 1);

for j = 1 : n - 1
    Q(:, j) = [ real(G{j}(1,1)) ; imag(G{j}(1,1)) ; real(G{j}(2,1)) ];
end

D = [ real(d).' ; imag(d).' ];

D = unitary_qr_eiscor(Q, D);
d = D(1, :) + 1i * D(2, :); d = d.';

its = 0;

end
