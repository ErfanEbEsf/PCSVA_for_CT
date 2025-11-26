function [dy] = dym(u)
% Implementation of finite differences
%each row - previous row, replicative boundary

[M,N] = size(u);
dy = [u(1:end-1,:,:);zeros(1,N)] - [zeros(1,N);u(1:end-1,:)];