fprintf('If necessary, I will download and install eiscor. \n');
fprintf('You need to have the following packages installed: git, gfortran, make. \n\n');
fprintf('Should I proceed? '); res = input('[yn]: ', 's');

if strcmp(res, 'y') || strcmp(res, 'Y')
    % Just in case unitary_qr_eiscor.mex* is here, we might need to clear
    % it now to avoid crashing MATLAB
    clear unitary_qr_eiscor

    if exist('eiscor', 'dir') ~= 7
        fprintf('Cannot find the eiscor folder, downloading it from Github ...');

        status = system('git clone -q git://github.com/eiscor/eiscor');

        if status ~= 0
            error('Failed to download eiscor! Check your git installation and internet connection');
        else
            fprintf(' done\n');
        end
    end

    if exist('~/eiscor/lib/libeiscor.so.0.2.0', 'file') ~= 2
      fprintf('EISCOR is apparently not installed, compiling and intalling into ~/eiscor');

      status = system('cd eiscor && make && make install && cd ..');

      if status ~= 0
            error('Failed to compile eiscor! Check your make and GNU Fortran installation');
        end
    end

    mex unitary_qr_eiscor.F90 ~/eiscor/lib/libeiscor.so.0.2.0
else
    fprintf('Ok, aborting.\n');
end
