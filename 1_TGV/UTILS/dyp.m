function [dy] = dyp(u)
% Implementation of finite differences
% each row - previos row, zero boundary

[M N P] = size(u);
dy = [u(2:end,:,:);u(end,:,:)] - u;
