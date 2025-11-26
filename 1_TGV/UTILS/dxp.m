function [dx] = dxp(u)
% Implementation of finite differences
% each col - previos col, zero boundary

dx = [u(:,2:end,:),u(:,end,:)] - u;
