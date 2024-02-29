%Muhammed Saadeddin Koçak 2232346

clearvars -except fanbeam sizeofimage D

projectionmatrix=fanbeam.projections;

[stepsize numberofbeams]=size(projectionmatrix);
imagematrix=zeros(sizeofimage);

%Loading Image
%Realize original image is loaded only for comparison. In backprojection
%code A matrix is never used.
myimage=load("square_circle.mat");
C=struct2cell(myimage);%Cell
A=cell2mat(C);%Matrix
figure
imagesc(A);
title ("Original Image");
colormap gray;

%Modifying Projection Matrix
for i=1:numberofbeams
    for j=1:stepsize
        modifiedprojectionmatrix(j,i)=projectionmatrix(j,i)*D*cosd((180/stepsize)*(j-(stepsize/2)));
    end
end
projectionmatrix=modifiedprojectionmatrix;

%Filtering
%Triangular Filter
triangfilter=triang(stepsize);%Designing filter in Fourier Domain

for o=1:numberofbeams
filteredprojectionmatrix(:,o)=ifft2((triangfilter).*(fft2(projectionmatrix(:,o))));
end
projectionmatrix=filteredprojectionmatrix;

%Defining beta, gamma, x, and y arrays according to beam number and step size
projectionanglerange=fanbeam.beta_range;
beta=linspace(projectionanglerange(1),projectionanglerange(2)-(projectionanglerange(2)/numberofbeams),numberofbeams);
gamma=linspace(-90,90-(180/stepsize),stepsize);


x=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);
y=linspace(-sizeofimage/2,sizeofimage/2,sizeofimage+1);

 for i=1:(length(beta))
     for j=1:(length(gamma))
%Clearing dynamic memory
relevantmatrix=[];
rowdata=[]; columndata=[]; midpointx=[];
midpointy=[]; distance=[];
%Finding intersection points
    yfromx=(D*sind(gamma(j))/sind(beta(i)+gamma(j)))-(x.*cotd(beta(i)+gamma(j)));
    xfromy=(D*sind(gamma(j))/cosd(beta(i)+gamma(j)))-(y.*tand(beta(i)+gamma(j)));
%Detecting relevant x and y coordinates
m=1;
for g=1:(length(x))
    if (((-sizeofimage/2)<=yfromx(g))&&(yfromx(g)<=sizeofimage/2))
        relevantmatrix(m,1)=x(1,g);
        relevantmatrix(m,2)=yfromx(1,g);
        m=m+1;
    else
    end
end
for k=1:(length(y))
    if (((-sizeofimage/2)<=xfromy(k))&&(xfromy(k)<=sizeofimage/2))
        relevantmatrix(m,1)=xfromy(1,k);
        relevantmatrix(m,2)=y(1,k);
        m=m+1;
    else
    end
end
uniquerelevantmatrix=unique(relevantmatrix,'rows');%Deleting same x and y points
sortedrelevantmatrix=sortrows(uniquerelevantmatrix);%Sorting relevant x and y points
sizerelevant=size(uniquerelevantmatrix,1);%Detecting if there are relevant points
if sizerelevant==0 %If no relevant points do nothing
else
[rowsizeofURM,columnsizeofURM]=size(uniquerelevantmatrix);
if (rowsizeofURM>1)%If there are more than one relevant point then calculate attenuation
%Calculating distance and midpoints to use them for detecting row and
%column data
for r=1:(rowsizeofURM-1)
    distance(r)=sqrt(((uniquerelevantmatrix(r,1)-uniquerelevantmatrix(r+1,1))^2)+((uniquerelevantmatrix(r,2)-uniquerelevantmatrix(r+1,2))^2));
    midpointx(r)=(uniquerelevantmatrix(r,1)+uniquerelevantmatrix(r+1,1))/2;
    midpointy(r)=(uniquerelevantmatrix(r,2)+uniquerelevantmatrix(r+1,2))/2;
end
%Finding corresponding column and row data
rowdata=(sizeofimage/2)-floor(midpointy);
columndata=(sizeofimage/2)+ceil(midpointx);
%Backprojecting one ray beam
sizeofdistancematrix=size(distance);
sizeofdistancearray=sizeofdistancematrix(1,2);
for h=1:sizeofdistancearray
    imagematrix(rowdata(1,h),columndata(1,h))=imagematrix(rowdata(1,h),columndata(1,h))+(projectionmatrix(j,i)*distance(1,h));

end
else%If there is only one relevant point beam is tangential, then do nothing
end
end
    end
 end
%Normalizing image matrix
imagematrix=real(imagematrix);
normimagematrix=(imagematrix-min(min(imagematrix)))/(max(max(imagematrix)));
%Displaying Image
figure
imagesc(normimagematrix);
title("Backprojected Image");
colormap gray;

%Quantitative Comparison
for i=1:sizeofimage
    for j=1:sizeofimage
        difference(i,j)=abs(normimagematrix(i,j)-A(i,j));
    end
end
%Displaying Differences Between Images
figure
imagesc(difference);
title("Difference Between Images");
colormap gray;
%Calculating Error
error=sum(sum(difference));
errorratio=error/(sizeofimage^2);