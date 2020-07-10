 function DKI_compute(optionFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   diff_kurt_tensor_parallel(optionFile)
%
%   ex: diff_kurt_tensor_parallel('DKIopt.mat');
% 
%   Computes the diffusion and kurtosis values for a data set. Results are
%   saved into the file name specified. Final results include lambda1,
%   lambda2, lambda3, vec1, vec2, vec3, K1, K2, K3, MK, exitLog, resLog, 
%   R2Log.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('E:\OneDrive\Frenzy_DKIoutput.mat');

readType = 'uint16';

numMatPool = options.numMatPool;
dim = options.dim;
b_value = options.b_value;
fname = options.fname;
gradientFile = options.gradientFile;
maskname = options.maskname;
S0thresh = options.S0thresh;
gradZeroIdx = options.gradZeroIdx;
outputFile = options.outputFile;
numB = numel(b_value);
vox = options.vox;
FWHM = options.FWHM;
maxK = 3;

%currMatPools = matlabpool('size');
%if currMatPools~=numMatPool
%    if currMatPools~=0
%        matlabpool close;
%    end
%    matlabpool(numMatPool);
%end
parpool(6)
fid = zeros(numB,1);
for bidx = 1:numB
    fid(bidx) = fopen(fname{bidx},'rb');    
end
fid
if isempty(maskname)
    mask = ones(dim(1),dim(2),dim(3));
else
    fidmask = fopen(maskname,'r');
    mask = fread(fidmask,Inf,'uint8');
    mask = reshape(mask,dim(1:3));    
end

timeLog = zeros(dim(3),1);

DWI = zeros([dim numB],'single');
lambda1 = zeros(dim(1),dim(2),dim(3),'single');
lambda2 = zeros(dim(1),dim(2),dim(3),'single');
lambda3 = zeros(dim(1),dim(2),dim(3),'single');

vec1 = zeros(3,dim(1),dim(2),dim(3),'single');
vec2 = zeros(3,dim(1),dim(2),dim(3),'single');
vec3 = zeros(3,dim(1),dim(2),dim(3),'single');

K1 = zeros(dim(1),dim(2),dim(3),'single');
K23 = zeros(dim(1),dim(2),dim(3),'single');
K23_simp = zeros(dim(1),dim(2),dim(3),'single');
MK = zeros(dim(1),dim(2),dim(3),'single');
ISO = zeros(dim(1),dim(2),dim(3),'single');

k23_temp = zeros(1,dim(1)*dim(2),'single');

b = zeros(numB*(dim(4)-1),21);
for bidx = 1:numB
    [btemp, Ad, Ak] = grad2bMx(gradientFile,b_value(bidx),gradZeroIdx);
    b((bidx-1)*(dim(4)-1)+1:bidx*(dim(4)-1),:) = btemp;
    if bidx==1
        bd = Ad;
    end
end

maxB = max(b_value);

gradIdxK = [1 1 1 1;
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
       
gradIdxD = [1 1;
            2 2;
            3 3;
            1 2;
            2 3;
            1 3];
            
dirComb = repmat([1 4 4 6 12 6 4 12 12 4 1 4 6 4 1]',1,dim(1)*dim(2));

Ad = -Ad/maxB;
Ak = Ak*6/maxB^2;

% Constraints matrix
maxC = 3;
C = [-Ad zeros(dim(4)-1,15); zeros(dim(4)-1,6) -Ak; -(maxC/maxB)*Ad Ak];

optimize = optimoptions('lsqlin');
optimize = optimoptions(optimize,'Algorithm','interior-point','Display','off','MaxIter',1500);


% Read in DWI's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bidx = 1:numB    
    a = fread(fid(bidx),Inf,[readType '=>single']);    
    a = reshape(a,dim);
    b0mask = a(:,:,:,gradZeroIdx)>S0thresh;
    DWI(:,:,:,:,bidx) = a;    
end
DWIorig = DWI;

% Gaussian Filter DWI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Spatially smoothing diffusion weighted images...');
f = genGaussKernel(vox,FWHM);
for bidx = 1:numB
    a = DWI(:,:,:,:,bidx);
    parfor g = 1:dim(4)
        a(:,:,:,g) = imfilter(a(:,:,:,g),f).*b0mask;
    end
    DWI(:,:,:,:,bidx) = a;
end
fclose all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Compute DKI metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
disp('Computing DKI metrics...');
sRange = (1:dim(3));
% sRange = 18;
for slice = sRange
    tic    
   
    % Display time elapsed and time remaining for current subject
    if slice~=sRange(1)
        disp(['Time for slice ' num2str(slice-1) ': ' num2str(timeLog(slice-1)) ' sec']);
        avgTimePerSlice = timeLog(slice-1);
        estRemTime = (sRange(end)-slice+1)*avgTimePerSlice;
        hrs = floor(estRemTime/3600);
        minutes = floor((estRemTime-hrs*3600)/60);
        sec = estRemTime-hrs*3600-minutes*60;
        disp(['Estimated remaining time: ' num2str(hrs) 'hrs   ' num2str(minutes) 'min   ' num2str(sec) 'sec']);
        elapsedTime = sum(timeLog(sRange(1):slice-1));
        hrs = floor(elapsedTime/3600);
        minutes = floor((elapsedTime-hrs*3600)/60);
        sec = elapsedTime-hrs*3600-minutes*60;
        disp(['Elapsed time: ' num2str(hrs) 'hrs   ' num2str(minutes) 'min   ' num2str(sec) 'sec']);
        disp(' ');
        disp('--------------------------------------------------------')
        disp(['        Number of Matlab Pools: ' num2str(numMatPool)]);
        disp(['               File Dimensions: ' num2str(dim(1)) ' X ' num2str(dim(2)) ' X ' num2str(dim(3)) ' X ' num2str(dim(4))]);
        for bidx = 1:numB
            eval(['disp([''               REC File Name ' num2str(bidx) ': '' fname{' num2str(bidx) '}]);']);
            eval(['disp([''                     b-value ' num2str(bidx) ': '' num2str(b_value(' num2str(bidx) '))]);']);
        end        
        disp(['            Gradient File Name: ' gradientFile]);        
        disp(['                Mask File Name: ' maskname]);
        disp(['0-Gradient Intensity Threshold: ' num2str(S0thresh)]);
        disp(['              0-Gradient Index: ' num2str(gradZeroIdx)]);
        disp(['              Output File Name: ' outputFile]);        
    end
    
    S = zeros(dim(1),dim(2),dim(4),numB);
    Sorig = zeros(dim(1),dim(2),dim(4),numB);
    S0 = zeros(dim(1),dim(2),numB);    
    sigRatio = cell(numB,1);
    tempIdx = ones(dim(4),1);
    tempIdx(gradZeroIdx) = 0;
    tempIdx = logical(tempIdx);
    for bidx = 1:numB                
        
        % mask DWI's for smoothed and original set
        for k = 1:dim(4)
            temp = DWI(:,:,slice,k,bidx);            
            S(:,:,k,bidx) = temp.*mask(:,:,slice);
            
            temp = DWIorig(:,:,slice,k,bidx);            
            Sorig(:,:,k,bidx) = temp.*mask(:,:,slice);
        end
        tempS0 = S(:,:,gradZeroIdx,bidx);
        tempS0(tempS0<S0thresh) = 0;
        
        % data set up for basic DTI computation on first b-value
        if (bidx==1)               
            ISO(:,:,slice) = mean(Sorig(:,:,tempIdx,bidx),3);
            tempS0orig = Sorig(:,:,gradZeroIdx,bidx);
            tempS0orig(tempS0orig<S0thresh) = 0;
            S0orig = tempS0orig;
            tempSbS0orig = log(Sorig(:,:,:,bidx)./repmat(S0orig,[1 1 dim(4)]));
            sigRatioMxorig = reshape(tempSbS0orig(:,:,tempIdx),dim(1)*dim(2),dim(4)-1); 
        end
        
        S0(:,:,bidx) = tempS0;        
        tempSbS0 = log(S(:,:,:,bidx)./repmat(S0(:,:,bidx),[1 1 dim(4)]));        
        sigRatio{bidx} = tempSbS0(:,:,tempIdx);                       
    end 
    
    sigRatioMx = zeros(dim(1)*dim(2),(dim(4)-1)*numB);
    for bidx = 1:numB        
        temp = sigRatio{bidx};
        sigRatioMx(:,(dim(4)-1)*(bidx-1)+1:(dim(4)-1)*bidx) = reshape(temp,dim(1)*dim(2),dim(4)-1);        
    end                
    
    % needed only for parfor due to indexing
    numGradDir = dim(4)-1;
    
    % loop through all voxels and calculate Diffusion and Kurtosis coefficients    
    compVoxIdx = zeros(dim(1)*dim(2),1);
    
    for pix = 1:dim(1)*dim(2)
        temp = sigRatioMx(pix,:);
        temp = temp(:);
        if sum(isinf(temp) | isnan(temp) | temp>1)==0 
            compVoxIdx(pix) = 1;
        end
    end
    disp(['#voxels to compute: ' num2str(sum(compVoxIdx))])
    compVoxIdx = logical(compVoxIdx);    
    tempSigRatioMx = sigRatioMx(compVoxIdx,:);
    X = zeros(21,dim(1)*dim(2));
    temp_X = zeros(21,sum(compVoxIdx));
    
    v = zeros(3,3,dim(1)*dim(2));
    d = zeros(3,3,dim(1)*dim(2));
    
    % compute DTI metrics using original method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor pix2 = 1:dim(1)*dim(2)    
        sR = sigRatioMxorig(pix2,:);
        sR = sR(:);
        if sum(isinf(sR) + isnan(sR))==0   
            D_temp = regress(sR,bd);
            D = zeros(3,3);
            D(1,1) = D_temp(1);
            D(2,2) = D_temp(2);
            D(3,3) = D_temp(3);

            D(1,2) = D_temp(4);
            D(2,1) = D_temp(4);

            D(3,2) = D_temp(5);
            D(2,3) = D_temp(5);

            D(3,1) = D_temp(6);
            D(1,3) = D_temp(6);

            [v(:,:,pix2),d(:,:,pix2)] = svd(D);
        end
    end  
    
    l1 = reshape(d(1,1,:,:),dim(1),dim(2));
    l2 = reshape(d(2,2,:,:),dim(1),dim(2));
    l3 = reshape(d(3,3,:,:),dim(1),dim(2));
    
    lambda1(:,:,slice) = l1;
    lambda2(:,:,slice) = l2;
    lambda3(:,:,slice) = l3;
  
    v_temp = reshape(v,9,dim(1)*dim(2));    
    for comp = 1:3
        vec1(comp,:,:,slice) = reshape(v_temp(comp,:),1,dim(1),dim(2));
        vec2(comp,:,:,slice) = reshape(v_temp(comp+3,:),1,dim(1),dim(2));
        vec3(comp,:,:,slice) = reshape(v_temp(comp+6,:),1,dim(1),dim(2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    parfor pix1 = 1:sum(compVoxIdx)      
        sigRatio = tempSigRatioMx(pix1,:);
        sigRatio = sigRatio(:);           
        % Constrained Linear Least Square Fitting (CLLS)
        [temp_X(:,pix1),~,~,~] = lsqlin(b,sigRatio,C,[zeros(numGradDir,1); zeros(numGradDir,1); zeros(numGradDir,1)],[],[],[],[],[],optimize);    
    end  
       
    X(:,compVoxIdx) = temp_X;
            
    % create diffusion tensor
    D = zeros(3,3,dim(1)*dim(2));
    D(1,1,:) = X(1,:);
    D(2,2,:) = X(2,:);
    D(3,3,:) = X(3,:);

    D(1,2,:) = X(4,:);
    D(2,1,:) = X(4,:);

    D(3,2,:) = X(5,:);
    D(2,3,:) = X(5,:);

    D(3,1,:) = X(6,:);
    D(1,3,:) = X(6,:);
        
    % calculate eigenvalues and eigenvectors of the diffusion tensor    
    parfor pix2=1:dim(1)*dim(2)
        [v(:,:,pix2),d(:,:,pix2)] = svd(D(:,:,pix2));
    end
    
    v_temp = reshape(v,9,dim(1)*dim(2));    
    v1 = v_temp(1:3,1:dim(1)*dim(2));
    v2 = v_temp(4:6,1:dim(1)*dim(2));
    v3 = v_temp(7:9,1:dim(1)*dim(2));
    
    Dapp = abs(Ad*X(1:6,:));
    Vapp = abs(Ak*X(7:21,:));
    Kapp = Vapp./(Dapp.^2);    
    
    Kapp(Kapp>maxK) = maxK;  
    Kapp(Kapp<0) = 0;

    MK(:,:,slice) = reshape(mean(Kapp),dim(1),dim(2));
    
    % find the AKC in the direction of the 3 eigenvectors
    Ak1111 = dirComb.*v1(gradIdxK(:,1),:).*v1(gradIdxK(:,2),:).*v1(gradIdxK(:,3),:).*v1(gradIdxK(:,4),:);
    Ak2222 = dirComb.*v2(gradIdxK(:,1),:).*v2(gradIdxK(:,2),:).*v2(gradIdxK(:,3),:).*v2(gradIdxK(:,4),:);
    Ak3333 = dirComb.*v3(gradIdxK(:,1),:).*v3(gradIdxK(:,2),:).*v3(gradIdxK(:,3),:).*v3(gradIdxK(:,4),:);

    % compute RK using by integrating around lambda 2,3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = atan2(v1(2,:),v1(1,:));
    phi = acos(v1(3,:));

    N = 360;
    dirCombK = repmat(dirComb(:,1),1,N);
    dirCombD = repmat([1 1 1 2 2 2]',1,N);
    
    for k = 1:dim(1)*dim(2)        
        ringVec = [cos(2*pi*(0:N-1)/N); sin(2*pi*(0:N-1)/N); zeros(1,N)];
        Ry = [cos(phi(k)) 0 sin(phi(k)); 0 1 0; -sin(phi(k)) 0 cos(phi(k))];
        Rz = [cos(theta(k)) -sin(theta(k)) 0; sin(theta(k)) cos(theta(k)) 0; 0 0 1];
        ringVec = Rz*Ry*ringVec;
        Akrad = dirCombK.*ringVec(gradIdxK(:,1),:).*ringVec(gradIdxK(:,2),:).*ringVec(gradIdxK(:,3),:).*ringVec(gradIdxK(:,4),:);
        Adrad = dirCombD.*ringVec(gradIdxD(:,1),:).*ringVec(gradIdxD(:,2),:);  
        Vapp = abs(Akrad'*X(7:21,k));
        Dapp = abs(Adrad'*X(1:6,k));
        k23_temp(k) = mean(Vapp./Dapp.^2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Wrot1111 = abs(sum(Ak1111.*X(7:21,:)));
    Wrot2222 = abs(sum(Ak2222.*X(7:21,:)));
    Wrot3333 = abs(sum(Ak3333.*X(7:21,:)));    
    d1 = reshape(d(1,1,:),1,dim(1)*dim(2));
    d2 = reshape(d(2,2,:),1,dim(1)*dim(2));
    d3 = reshape(d(3,3,:),1,dim(1)*dim(2));
    k1_temp = Wrot1111./d1.^2;
    k2_temp = Wrot2222./d2.^2;
    k3_temp = Wrot3333./d3.^2;
    
    k1 = max(min(k1_temp,maxK),0);
    k23 = max(min(k23_temp,maxK),0);
    k23_simp = max(min((k2_temp+k3_temp)/2,maxK),0);
    
    k1(~compVoxIdx) = 0;
    k23(~compVoxIdx) = 0;
    k23_simp(~compVoxIdx) = 0;
    
    K1(:,:,slice) = reshape(k1,dim(1),dim(2));
    K23(:,:,slice) = reshape(k23,dim(1),dim(2));
    K23_simp(:,:,slice) = reshape(k23_simp,dim(1),dim(2));
      
    
    timeLog(slice) = toc;
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
end
 
MK(isnan(MK)) = 0;
K1(isnan(K1)) = 0;
K23(isnan(K23)) = 0;
K23_simp(isnan(K23_simp)) = 0;

fclose all;

rawFilename = regexp(outputFile,'\.','split');
rawFilename = rawFilename{1};

fidKmean = fopen([rawFilename '.MK'],'wb');
fidKaxial = fopen([rawFilename '.AK'],'wb');
fidKradial = fopen([rawFilename '.RK'],'wb');
fidKradial2 = fopen([rawFilename '.RK2'],'wb');
fidDmean = fopen([rawFilename '.MD'],'wb');
fidDaxial = fopen([rawFilename '.AD'],'wb');
fidDradial = fopen([rawFilename '.RD'],'wb');
fidFA = fopen([rawFilename '.FA'],'wb');
fidVEC1 = fopen([rawFilename '.V1'],'wb');
fidISO = fopen([rawFilename '.ISO'],'wb');

MD = (lambda1+lambda2+lambda3)/3;
FA = cal_fa(lambda1,lambda2,lambda3);

fwrite(fidFA,FA,'single');
fwrite(fidKmean,MK,'single');
fwrite(fidKaxial,K1,'single');
fwrite(fidKradial,K23,'single');
fwrite(fidKradial2,K23_simp,'single');
fwrite(fidDmean,MD,'single');
fwrite(fidDaxial,lambda1,'single');
fwrite(fidDradial,(lambda2+lambda3)/2,'single');
fwrite(fidVEC1,vec1,'single');
fwrite(fidISO,ISO,'single');
fclose all;

save(outputFile,'lambda1','lambda2','lambda3','vec1','vec2','vec3','K1','K23','K23_simp','MK','ISO');
