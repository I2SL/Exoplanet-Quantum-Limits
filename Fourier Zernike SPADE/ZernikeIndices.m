function [nj,mj] = ZernikeIndices(n_max)
    % returns the radial and azimuthal indices up to radial
    % order n_max using the 
    nj = [];
    mj = [];
    for n = 0:n_max
        for m = -n:2:n
            nj = [nj, n];
            mj = [mj, m];
        end
    end
end