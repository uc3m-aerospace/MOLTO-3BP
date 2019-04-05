function Mat = VectoMat6(Vec)

% This function converts a 1x36 vector (Vec) in a 6x6 matrix (Mat)

Mat=zeros(6);
for i=1:6
    for j=1:6
        k=j+6*(i-1);
        Mat(i,j)=Vec(k);
    end
end