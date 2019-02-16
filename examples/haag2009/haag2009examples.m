% examples from
% Thomas Haag, Ulrich Münz, Frank Allgöwer:
% Comparison of Different Stability Conditions for Linear Time-Delay
% Systems with Incommensurate Delays, 2009

% dx = A0*x(t) + \sum_{k=1}^N A_k*x(t-\tau_k)
%
% problem: find maximum delay \bar{\tau} so that all \tau_k < \bar{\tau}

clear all;

filename_base = 'haag2009ex';
filename_ext  = 'chEq';

% Example 1

A0 = [0 1; -1 -1];
A1 = [0 0; 0 -1];
c1 = characteristic_eq(A0,A1);

% Example 2

A0 = [-2 0; 0 -0.9];
A1 = [-1 0; -1 -1];
c2 = characteristic_eq(A0,A1);

% Example 3

A0 = [-1 13.5 -1; -3 -1 -2; -2 -1 -4];
A1 = [-5.9 7.1 -70.3; 2 -1 5; 2 0 6];
c3 = characteristic_eq(A0,A1);

% Example 4

A0 = [0.5];
A1 = [-0.9];
A2 = [-1.5];
c4 = characteristic_eq(A0,A1,A2);

% Example 5

A0 = [0];
A1 = [-1];
A2 = [-2];
c5 = characteristic_eq(A0,A1,A2);

% Example 6

A0 = [0 1; -(1+(49/256)) -7/8];
A1 = [0 0; 1/5 0];
A2 = [0 0; -4/5 0];
c6 = characteristic_eq(A0,A1,A2);

% Example 7

A0 = [-1 13.5 -1; -3 -1 -2; -2 -1 -4];
A1 = [-5.9 0 0; 2 0 0; 2 0 0];
A2 = [0 7.1 -70.3; 0 -1 5; 0 0 6];
c7 = characteristic_eq(A0,A1,A2);

% Example 8

A0 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -2 -3 -5 -2];
A1 = [-9/200 1.5/200 1/4 0; 1/200 1/200 1/20 0; 0 0 0 1/2000; -2 -1/2 -1 0];
A2 = [3/80 0 3/40 1/8; 0 1/20 1/20 0; 1/20 1/20 0 0; 0 -2.5 0 -1];
c8 = characteristic_eq(A0,A1,A2);

% Example 9

A0 = [-2];
A1 = [-4];
A2 = [-3];
A3 = [-1];
c9 = characteristic_eq(A0,A1,A2,A3);

% store to files
for i=1:9
    fid = fopen([filename_base num2str(i) '.' filename_ext], 'w');
    try
    eval(['fprintf(fid, ''%s'', [c' num2str(i) ']);']);
    catch err
        warning(['Could not write c' num2str(i) ' to file!']);
    end
    fclose(fid);
end