function tm = fast_timeit(f)
%FAST_TIMEIT 

% We try tic / toc first, and resort to timeit() only if the time is small
% enough. 
timer = tic;
f();
tm = toc(timer);

if tm < 1
    tm = timeit(f);
end

end

