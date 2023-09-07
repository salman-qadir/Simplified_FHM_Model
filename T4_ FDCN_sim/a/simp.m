function I = simp(f)
% This function computes the integral "I" via Simpson's rule  with equally spaced points
% Syntax: I = simpsons(f) 
% Where,
%  f= can be either an anonymous function (e.g. f=@(x) sin(x)) or a vector
%  containing equally spaced values of the function to be integrated
% Example
% Suppose you want to integrate a function f(x) in the interval [-1,1].
% You know some values of the function f(x) between the given interval,
% those are fi= {1,0.518,0.230,0.078,0.014,0,0.006,0.014,0.014,0.006,0}
% Thus:
% fi= [1 0.518 0.230 0.078 0.014 0 0.006 0.014 0.014 0.006 0];
% I=simpsons(fi)
if numel(f)>1 % If the input provided is a vector

    I= 1/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end))+f(end));
end
