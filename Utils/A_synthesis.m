function [ y_vector ] = A_synthesis( u_vector, w, S)
%A_synthesis takes in VECTOR of floor coefficients, u, and constructs the floor image
%using the basis specified by w, to the level specified by l
%u is the image RESHAPED into a rectangle!

y = waverec2(u_vector, S, w);
% reshape y into a column
y_vector = y(:);

end

