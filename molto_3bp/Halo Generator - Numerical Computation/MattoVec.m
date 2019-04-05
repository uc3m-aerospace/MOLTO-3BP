function Vec = MattoVec(Mat)

% This function converts a nxn matrix (Mat) in a 1xn^2 vector (Vec)
    % The vector from 1 to n^2 is formed by the rows of the matrix

[r,c]=size(Mat);
Vec=zeros(r*c,1);
for i=1:r
    for j=1:c
        k=j+r*(i-1);
        Vec(k)=Mat(i,j);
    end
end