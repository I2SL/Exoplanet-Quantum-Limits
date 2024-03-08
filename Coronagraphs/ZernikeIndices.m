function [nj,mj] = ZernikeIndices(n_max)
    % returns the radial and azimuthal indices up to radial
    % order n_max using the OSA/ANSI indexing standard 
    % (up to shift by +1 for matlab indexing)
    %
    % Author(s): Nico Deshler, University of Arizona
    % Affiliation(s): Wyant College of Optical Sciences, University of Arizona
    % Date: March 7, 2024
    
    nj = [];
    mj = [];
    for n = 0:n_max
        for m = -n:2:n
            nj = [nj, n];
            mj = [mj, m];
        end
    end
end