% examples from
% Fabien Lescher and Clement Roos:
% Robust stability of time-delay systems with structured uncertainties: a mu-analysis based algorithm
% CDC-ECC, 2011

% dx = A0*x(t) + \sum_{k=1}^N A_k*x(t-\tau_k)
%
% problem: find maximum delay \tau so that system is stable for all
% delta or all delta_i=[-1,1]

clear all;

filename_base = 'lescher2011ex';
filename_ext  = 'chEq';

% Example 1
syms delta real

A0 = [0 -0.12+0.42*delta; 1 -0.465-0.035*delta];
A1 = [-0.1 -0.35; 0 0.3];
c1 = characteristic_eq(A0,A1);

% Example 2 (adapted from haag2009 or references within, was ex2 there as well but without deltas!)
syms delta_1 delta_2 delta_3 delta_4 real

A0 = [-2+1.6*delta_1 0; 0 -0.9+0.05*delta_2];
A1 = [-1+0.1*delta_3 0; -1 -1+0.3*delta_4];
c2 = characteristic_eq(A0,A1);

% store to files
for i=1:2
    fid = fopen([filename_base num2str(i) '.' filename_ext], 'w');
    try
    eval(['fprintf(fid, ''%s'', [c' num2str(i) ']);']);
    catch err
        warning(['Could not write c' num2str(i) ' to file!']);
    end
    fclose(fid);
end