function [ fin ] = findc(p,x)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fin=find(min(abs(p-x))==abs(p-x));

end

