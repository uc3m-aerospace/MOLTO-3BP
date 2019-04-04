function [x0, y0, vx0, vy0, eigvec, eigval, inv_phi_0] = IC_state_matrix(Ax,a,b)

% Compute state transition matrix
A = [0 0 1 0; 0 0 0 1; a 0 0 2; 0 -b -2 0];
[eigvec,eigval] = eig(A);
eigval = diag(eigval);%returns a column vector with the main diagonal 
Mat = [eigvec(1,3)+eigvec(1,4)];
inv_phi_0 = inv(eigvec);
x0 = Ax;%Initial position in the x axis at a distance x0 from Libration point
y0 = 0;
initpos = x0;

coef = Mat\initpos;

c1 = 0;
c2 = 0;
c3 = coef;
c4 = c3;

vx0 = c1*eigvec(3,1)+c2*eigvec(3,2)+c3*eigvec(3,3)+c4*eigvec(3,4);%=0
vy0 = c1*eigvec(4,1)+c2*eigvec(4,2)+c3*eigvec(4,3)+c4*eigvec(4,4);
end