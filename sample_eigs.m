function e = sample_eigs(n, complexoutput, detsign, witheiscor)
%SAMPLE_EIGS Sample the eigenvalues of Haar distributed matrices.
%
%   E = SAMPLE_EIGS(N,COMPLEXOUTPUT) are the eigenvalues sampled from the Haar
%   distribution over the unitary group U(N), if COMPLEXOUTPUT is TRUE, or over
%   the orthogonal group O(N), if COMPLEXOUTPUT is FALSE. If the second argument
%   is not provided, COMPLEXOUTPUT defaults to TRUE.
%
%   E = SAMPLE_EIGS(N,COMPLEXOUTPUT,DETSIGN) sample from matrices with a
%   specified determinant. If COMPLEXOUTPUT is TRUE, then the algorithms samples
%   matrices with determinant EXP(ANGLE(DETSIGN)*I)), whereas if COMPLEXOUTPUT
%   is FALSE the algorithms samples matrices with determinant SIGN(DETSIGN). An
%   empty DETSIGN is ignored, and an error is produced if DETSIGN is 0.
%   Setting DETSING to 1 samples from the Haar distribution over the special
%   unitary group SU(n), if COMPLEXOUTPUT is TRUE, and from the Haar
%   distribution over the supecial orthogonal group SO(n) if COMPEXOUTPUT is
%   FALSE.
%
%   E = SAMPLE_EIGS(N,COMPLEXOUTPUT,DETSIGN,WITHEISCOR) uses the eiscor package
%   (https://github.com/eiscor/eiscor) if the eiscor library if WITHEISCOR is
%   TRUE and the mex file unitary-qr-eiscor is in the MATLAB path, otherwie it
%   uses a pure MATLAB implementation. Deafult value is TRUE.
%
%   Reference:
%
%   [1] M. Fasi & L. Robol. Sampling the eigenvalues of random orthogonal and
%   unitary matrices. In preparation.

% Parse input.
  if (nargin < 1 || nargin > 4)
    error('sampleeigs:invalidArgumentNumber',...
          'This function accepts only one to four arguments.');
  end
  if ~(isnumeric(n) && isscalar(n) && round(n) == n && n > 0)
    error('sampleeigs:invalidN',...
          'The first argument must be a positive scalar integer.');
  end
  if nargin < 2
    complexoutput = true;
  else
    if ~(islogical(complexoutput) && isscalar(complexoutput))
      error('sampleeigs:invalidComplexoutput',...
            'The second argument must be logical.');
    end
  end
  fixsign = false;
  if nargin >= 3 && ~isempty(detsign)
    if ~(isnumeric(detsign) && isscalar(detsign)) || detsign == 0
      error('sampleeigs:invalidN',...
            'The third argument must be a non-zero scalar.');
    else
      fixsign = true;
    end
  end
  if nargin < 4
    witheiscor = true;
    
    % Check that eiscor actually exists
    if exist('unitary_qr_eiscor', 'file') ~= 3
        warning('sampleeigs:noeiscor', ...
            'EISCOR is not available: relying on the standard QR solver of MATLAB');
        witheiscor = false;
    end
  else
    if ~(islogical(witheiscor) && isscalar(witheiscor))
      error('sampleeigs:invalidComplexoutput',...
            'The second argument must be logical.');
    end
  end

  % Generate Hessenberg form of Haar distributed matrix.
  u = [randn(1,n-1);
       chi2rnd(n-1:-1:1)];
  if complexoutput
    u(1,:) = u(1,:) + 1i*randn(1,n-1);
    u(2,:) = u(2,:) + chi2rnd(n-1:-1:1); % = chi2rnd(2*(n-1:-1:1));
  end
  u(2,:) = sqrt(u(2,:));
  if complexoutput
    dn = randn(1) + 1i*randn(1);
  else
    dn = randn(1);
  end
  d = -[sign([u(1,:), dn])];
  norms = sqrt(sum(abs(u).^2,1));
  u(1,:) = u(1,:) - d(1:n-1) .* norms;

  % D = diag(d);
  % for j = n-1:-1:1
  %   U = eye(2) - 2/(norm(u(:,j))^2) * u(:,j) * u(:,j)';
  %   D(j:j+1,:) = U * D(j:j+1,:);
  % end
  % e = eig(D);

  % return

  % Refactor core transformations.
  Q = zeros(3,n-1);
  d2 = 1;
  factors = -2./sum(abs(u).^2,1);
  for j = 1:n-1
    % Build 2-by-2 Householder matrix.
    % U = eye(2) - 2/norm(u(:,j))^2 * u(:,j) * u(:,j)';
    U = (factors(j) * u(:,j)) * u(:,j)';
    U([1,4]) = U([1,4]) + 1;
    % Scale first row.
    U(1,:) = d2*U(1,:);
    sgn = sign(U(2,1)');
    % constant = 2./(U(1,1).^2+U(2,2).^2);
    c(j) = sgn * U(1,1);
    s(j) = -sgn*U(2,1);
    d(j) = d(j) * (conj(sgn)*abs(U(1,1))^2+sgn*U(2,1)^2);
    d2 = sgn * (U(1,1)*U(2,2)-U(2,1)*U(1,2)); % = sgn * det(U)
  end
  Q = [real(c); imag(c); real(s)];

  if fixsign
    if ~complexoutput
      d(n) = sign(real(detsign)) / prod(d(1:n-1));
    else
      d(n) = sign(detsign) / prod(d(1:n-1));
    end
  else
    d(n) = d(n) * d2;
  end

  % Compute eigenvalues.
  if witheiscor
    % Make sure the diagonal is of unit elements up to machine precision,
    % eiscor can be quite picky about this fact. The same holds for Q
    d = d ./ abs(d);

    % Qmod = arrayfun(@(j) norm(Q(:,j)), 1 : size(Q, 2));
    Qmod = sqrt(sum(abs(Q).^2,1));
    Q = Q ./ ( ones(3, 1) * Qmod );

    e = unitary_qr_eiscor(Q, [ real(d) ; imag(d) ]);
    e = e(1, :) + 1i * e(2, :);
    e = e.';
  else
    D = diag(d);
    for j = n-1:-1:1
      D(j:j+1,:) = [c(j) s(j); -s(j) c(j)'] * D(j:j+1,:);
    end
    e = eig(D);
  end
end
