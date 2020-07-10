function [b, bd, bk] = grad2bMx(fileName, b_value, gradZeroIdx)

grad = csvread(fileName);
grad(:,3) = -grad(:,3);
gradValid = (1:size(grad,1))~=gradZeroIdx;
Gx = grad(gradValid,1);
Gy = grad(gradValid,2);
Gz = grad(gradValid,3);

bd = -b_value*[grad(gradValid,:).^2 2*Gx.*Gy 2*Gy.*Gz 2*Gx.*Gz];

gradIdx = [1 1 1 1;
           1 1 1 2;
           1 1 1 3;
           1 1 2 2;
           1 1 2 3;
           1 1 3 3;
           1 2 2 2;
           1 2 2 3;
           1 2 3 3;
           1 3 3 3;
           2 2 2 2;
           2 2 2 3;
           2 2 3 3;
           2 3 3 3;
           3 3 3 3];
       
dirComb = repmat([1 4 4 6 12 6 4 12 12 4 1 4 6 4 1],size(grad,1)-1,1);
       
bk = dirComb.*b_value^2/6.*grad(gradValid,gradIdx(:,1)).*grad(gradValid,gradIdx(:,2)).*grad(gradValid,gradIdx(:,3)).*grad(gradValid,gradIdx(:,4));

b = [bd bk];