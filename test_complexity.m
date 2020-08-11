% Test complexity of sample_eigs(n), wrt to n.

rng(0);
k = 15;

N = 2.^(1 : k);
tm = zeros(1, k);
tm_eig = zeros(1, k);

for j = 1 : k
    % Estimate the time needed to solve the problem. When n is small
    % enough, we use timeit() which is more accurate. If n > 4000, then we
    % resort to tic / toc which is much faster.

    n = N(j);

    tm(j) = fast_timeit(@() sample_eigs(n, true, [], true));

    % We avoid running eig on very large matrices
    if n <= 8192
      tm_eig(j) = fast_timeit(@() sample_eigs(n, true, [], false));
    else
      tm_eig(j) = NaN;
    end

    fprintf('N  = %d, time required = %f s, time for eig = %f s\n', n, tm(j), tm_eig(j));
end

loglog(N, tm, 'ro-', N, tm_eig, 'bx-');
legend('sample\_eigs', 'eig');
title('Complexity of sample\_eigs and standard eig');
xlabel('N');
ylabel('Time (s)');

dlmwrite('complexity.dat', [ N.', tm.' , tm_eig.'], '\t');