% test_eigenvalue_ditribution_spacing.m:
% Plot the distribution of the phases and normalized spacing of the
% eigenvalues of random unitary matrices sampled in different ways from the
% Haar distribution.

%% Initialization.
rng(193627)

nreps = 10000;
n = 10;

fprintf('TEST_DISTRIBUTION: Using nreps = %d and n = %d\n', nreps, n);
fprintf('                   Parameters can be adjusted editing test_distribution.m\n');

saveoutput = true;
useeiscor  = true;
useassert  = false;

% 4 configurations:
%  1) Complex unitary matrices
%  2) Orthogonal matrices
%  3) Matrices in O+
%  4) Matrices in O-

%% Compare eigenvalue distribution and spacing.
for j = 0 : 4
  switch j
    case 0
      usecomplex = true;
      dd = 1;
    case 1
      usecomplex = true;
      dd = [];
    case 2
      usecomplex = false;
      dd = [];
    case 3
      usecomplex = false;
      dd = 1;
    case 4
      usecomplex = false;
      dd = -1;
  end

  eigenvalues = zeros(1, n*nreps);
  distances = zeros(1, n*nreps);
  index_eigenvalues = 1;
  index_distances = 1;
  for i = 1:nreps

    e = sample_eigs(n, usecomplex, dd, useeiscor);

    % Compute phases of the eigenvalues.
    e = angle(e);
    e(e < 0) = e(e < 0) + 2*pi;
    if useassert
      assert(all(e >= 0))
      assert(all(e < 2*pi))
    end
    e = sort(e);
    eigenvalues(index_eigenvalues:index_eigenvalues+n-1) = e;

    % Compute spacing of the eigenvalues.
    e(n+1) = e(1)+2*pi;
    distances(index_distances:index_distances+n-1) = e(2:end) - e(1:end-1);
    index_eigenvalues = index_eigenvalues + n;
    index_distances = index_distances + n;
  end
  distances = n/(2*pi) * distances; % Normalization of the spacing.

  figure(j+1)
  clf
  
  subplot(2,1,1)
  hold on
  h1 = histogram(eigenvalues, 'normalization', 'pdf','BinLimits',[0 2*pi]);
  plot([0,2*pi],1/(2*pi) * ones(2,1),'r--','Linewidth', 2);
  hold off

  subplot(2,1,2)
  hold on
  h2 = histogram(distances, 'normalization', 'pdf');
  s = [0:0.01:3.5];
  plot(s, 32*(s/pi).^2 .* exp(-4*s.^2/pi), 'r--','Linewidth',2);
  hold off
  
  switch j
      case 0
          sgtitle('SU(n) - Complex unitary matrices with det = 1');
      case 1
          sgtitle('U(n) - Complex unitary matrices');
      case 2
          sgtitle('O(n) - Real orthogonal matrices');
      case 3
          sgtitle('SO(n) - Orthogonal matrices with det = 1');
      case 4
          sgtitle('O^-(n) - Orthogonal matrice with det = -1');
  end  

  if saveoutput
    if isempty(dd)
      dd = 0;
    end
    if usecomplex
      suffix = sprintf('unitary-%2d-%d', n, dd);
    else
      suffix = sprintf('orthog-%2d-%d', n, dd);
    end

    filename = sprintf('eig-dist-%s.dat', suffix);
    fid = fopen(filename, 'w+');
    data = [h1.BinEdges', [h1.Values'; h1.Values(end)]];
    for i = 1:length(h1.BinEdges)
      fprintf(fid, '%.15f %.15f\n', data(i,:));
    end
    fclose(fid);

    filename = sprintf('eig-spacing-%s.dat', suffix);
    fid = fopen(filename, 'w+');
    data = [h2.BinEdges', [h2.Values'; h2.Values(end)]];
    for i = 1:length(h2.BinEdges)
      fprintf(fid, '%.15f %.15f\n', data(i,:));
    end
    fclose(fid);
  end

end