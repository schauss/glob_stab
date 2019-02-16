function chEq = characteristic_eq(varargin)
% create characteristic equation for glob_stab from system matrices Ai
%
% give the matrices in the order A0, A1, A2, ...
% the system equation is dx(t) = A0*x(t) + A1*x(t-T_d1) + A2*x(t-T_d2) + ...

if length(varargin) < 1
    error ('No Input!');
end

A = varargin{1};

for i=2:length(varargin)
    if (size(A) ~= size(varargin{i}))
        error ('System matrices have different size!');
    end
    eval(['syms z' num2str(i-1) ';']);
    eval([' A = A + varargin{' num2str(i) '} * z' num2str(i-1) ';']);
end

%epsilon = 1e-10;
syms epsilon real
syms s;

p = det(eye(size(A))*(s-epsilon)-A);
p = collect(p,'s');
c = coeffs(p,'s');

chEq = [];

for i=2:length(varargin)
%     if isa(epsilon,'sym')
%         chEq = [chEq 'z' num2str(i-1) ' = exp(-T_d' num2str(i-1) '*(s-epsilon) );'];
%     elseif (epsilon == 0)
%         chEq = [chEq 'z' num2str(i-1) ' = exp(-T_d' num2str(i-1) '*s);'];
%     else
%         chEq = [chEq 'z' num2str(i-1) ' = exp(-T_d' num2str(i-1) '*(s-' num2str(epsilon) ') );'];
%     end
    chEq = [chEq 'z' num2str(i-1) ' = exp(-T_d' num2str(i-1) '*s)'];
    if isa(epsilon,'sym')
        chEq = [chEq '*exp(T_d' num2str(i-1) '*epsilon)'];
    elseif (epsilon ~= 0)
        chEq = [chEq '*exp(T_d' num2str(i-1) '*' num2str(epsilon) ')'];
    end
    chEq = sprintf('%s;\n',chEq);
end

cheq_t = ccode(c);
for i=1:length(c)
    cheq_t = strrep(cheq_t,['  t' num2str(i-1) ' ='],['chEq(' num2str(i) ') =']);
end
chEq = [chEq cheq_t];
chEq = sprintf('%s\n',chEq);