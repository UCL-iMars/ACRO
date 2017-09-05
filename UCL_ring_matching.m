function [RP, TP] = ring_matching(Sift_R_Descr, Sift_R_Coord, Sift_T_Descr, Sift_T_Coord, InpParam)

%Title: Ring Matching Matlab Software
%Developed by Dr Panagiotis Sidiropoulos, MSSL/UCL
%Version 1.0, 4-Sep-2017
%Contact: p.sidiropoulos@ucl.ac.uk
%
%External Dependencies: Andrea Vedaldi SIFT implementation: http://vision.ucla.edu/~vedaldi/code/sift.html
%
%
%Cite: P. Sidiropoulos and J.-P. Muller, A systematic solution to multi-instrument co-registration of high-resolution planetary images to an orthorectified baseline,
%      IEEE Transactions on Geoscience and Remote Sensing, 2017.
%
%
%Software Description: Matches points in two corresponding images using geometric constraints acording to the cited paper
%
%Software Input: SIFT-points and their coordinates of the reference (Sift_R_Descr, Sift_R_Coord) and the target (Sift_T_Descr, Sift_T_Coord) image plus a parameter
%vector. If the parameter vector is not defined default values are used. Both the SIFT-points and the coordinates are given as NXD matrices, where N is the number 
%of SIFT-points and D the dimension (D=128 for the SIFT-points and D=5 for their coordinates.
%
%Note: The coordinates file is a NX5 matrix, with the 5 columns being: the 2 first are the pixel-coordinates (point coordinates in pixels), the 3rd one is a 
%placeholder for future development and is currently irrelevant (default value is 0) and the 4-5th are the metre coordinates (point coordinates in metres, 
%according to some datum. Estimate the metre coordinates require georeferencing, but georeferecing is not relevant to the software (i.e. it is completed prior to it).
%
%Software Parameters:
%1: TimeThresh: Time (in seconds) that the software will continue looking for matches before aborting, default value 3600
%2: Ring_Radius: The radius of each ring (in kilometres), default value 4
%3: Max_Radius: The maximum distance (in kilometres) that each point can be matched, default value 40
%4: Target_Points: The number of points required for a ring to be identified as the correct on, default value 15
%5: Tolerance_Value: The tolerance in the outlier detection method (see the cited paper for more info), default value 0.02
%
%Software Output: The coordinates of matching points (in pixel coordinates) in the reference (RP) and the target (TP) image, respectively
%


if (nargin<4)||(nargin>5)
    error('Wrond number of input parameters');
elseif nargin==4
    InpParam = [3600 4 40 15 0.02];  
end
TimeThresh = InpParam(1);
Ring_Radius = InpParam(2);
Max_Radius = InpParam(3);
Target_Points = InpParam(4);
Tolerance_Val = InpParam(5);

NumT = size(Sift_T_Coord,1);
NumR = size(Sift_R_Coord,1);

if (NumT<Target_Points)||(NumR<Target_Points)
    error('Not Enough SIFT Points for the matching to proceed.');
elseif (size(Sift_R_Descr,2)~=128)||(size(Sift_T_Descr,2)~=128)||(size(Sift_R_Coord,2)~=5)||(size(Sift_T_Coord,2)~=5)
    error('Wrong format of SIFT points.');
end

F = 0; RP = []; TP = [];

Rad_Num = round(Max_Radius/Ring_Radius);
PointsNum = zeros(1,Rad_Num);
P = zeros(Target_Points,4*Rad_Num);

tic;
metr = 1 ;
while F==0
	%-------Break Check
    if NumR==1
        error('Run out of points, not enough matches found.');
    end
    
    T = toc;
    if T>TimeThresh
        error('Run out of time, not enough matches found.');
    end
    %-------Break Check
    
    %--------Find Matches
	R = randperm(NumR,1);
    Cur_R_Descr = Sift_R_Descr(R,:);
    Cur_R_Coord = Sift_R_Coord(R,1:2);
    Cur_R_Metr = Sift_R_Coord(R,4:5);
    
    Dist_Matrix = [Sift_T_Coord(:,4)-Cur_R_Metr(1) Sift_T_Coord(:,5)-Cur_R_Metr(2)];
    D = (Dist_Matrix(:,1).*Dist_Matrix(:,1)+Dist_Matrix(:,2).*Dist_Matrix(:,2)).^(0.5);
    
    MaxZ = 0;
    for i = 1 : Rad_Num
        Z = (D<(i*Ring_Radius)).*(D>((i-1)*Ring_Radius));
        MaxZ = max([MaxZ nnz(Z)]);
        if nnz(Z)>Target_Points
            Cur_T_Descr = Sift_T_Descr(Z==1,:);
            Cur_T_Coord = Sift_T_Coord(Z==1,1:2);
        
            cd sift;
            CurMatch = siftmatch(Cur_R_Descr', Cur_T_Descr', 1.5);
            cd ..
        
            if ~isempty(CurMatch)
                PointsNum(i) = PointsNum(i)+1;
                P(PointsNum(i),((4*i-3):4*i)) = [Cur_R_Coord Cur_T_Coord(CurMatch(2),:)];
            end
        end       
    end
    %--------Find Matches
    
    %--------Check Rings
    if mod(metr,200)==0
        Max_Cons = 3;
        Cand_Rings = (PointsNum>=Target_Points);
        for j = 1 : Rad_Num
            if Cand_Rings(j)
                X = ransac_permute(P(1:PointsNum(j),(4*j-3):(4*j-2)), P(1:PointsNum(j),(4*j-1):(4*j)), Tolerance_Val, Target_Points);
                Max_Cons = max([Max_Cons size(X,1)]);
                if size(X,1)>=Target_Points
                    F = 1;
                    RP = X(:,1:2);
                    TP = X(:,3:4);
                end
            end
        end
        %Max_Cons
    end
    %--------Check Rings
    
    Sift_R_Descr(R,:) = [];
    Sift_R_Coord(R,:) = [];
    NumR = size(Sift_R_Coord,1);
    metr = metr+1;
end

return



%---Ring matching outlier detection function
function Cur_P = ransac_permute(Ref_Control_Points, Tar_Control_Points, Tolerance_Val, Target_Points)

N = size(Ref_Control_Points,1);
Dist_Mat = zeros(N,N);
for i = 2 : N
    vr1 = Ref_Control_Points(i,:);
    vt1 = Tar_Control_Points(i,:);
    for j = 1 : (i-1)
        vr2 = Ref_Control_Points(j,:);
        vt2 = Tar_Control_Points(j,:);
        
        dr = vr2-vr1;
        dt = vt2-vt1;
        
        Dist_Mat(i,j) = sqrt(dr(1)^2+dr(2)^2)/sqrt(dt(1)^2+dt(2)^2);
        Dist_Mat(j,i) = Dist_Mat(i,j);
    end
end

Min_Lim = 1 - Tolerance_Val;
Max_Lim = 1 + Tolerance_Val;

Q1 = (Dist_Mat<Max_Lim);
Q2 = (Dist_Mat>Min_Lim);
Q = Q1.*Q2;

S = sum(Q,2);
if max(S)>Target_Points
    M = (Target_Points+1)*eye(N);
    for i = 2 : N
        for j = 1 : (i-1)
            M(i,j) = Q(i,:)*(Q(j,:)');
            M(j,i) = M(i,j);
        end
    end

    QM = (M>Target_Points);
    SQM = sum(QM);
    [~, F] = max(SQM);
    F2 = find(QM(F,:));
    Cur_P = [Ref_Control_Points(F2,:) Tar_Control_Points(F2,:)];
else
    Cur_P = [];
end

return
