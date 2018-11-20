%CaOx images 20x - segmentation optimization - noise reduction 2018/05/28 and manual classificaiton
%preparation of training data for CaOx classification

%create Training data with shape, intensity and texture features

%segmentation advice:
%reduce structuring element for wathershed based on intensities from 6
%-->5 for high particle number/touching images



%% ---INITIALIZATION---
clear variables
close all
clc

%% ---read original image---
cd 'D:\lab\Experiments\20180329_CaOxgrowthScreening-microscope\Matlab_Classification_protocols-evaluation\20180528_CaOxmicroscopy_classification09'%direct to folder
I = imread ('20180523_CaOxmiccroscopy_Exp06_ctrl_NaOx_w1_01_ch00.tif');%select input image
figure(1); clf
imshow (I);


%% ---detect edges---
Ismooth = imgaussfilt(I,1.4);
figure (1), imshow (Ismooth);
threshold = 0.03; % edge detection 
fudgeFactor = 0.5;
BWs = edge(Ismooth,'sobel', threshold * fudgeFactor); %make binary image of edges
figure (2), imshow(BWs), title('binary gradient mask');


%% --- morphologically close image and fill holes--
se = strel('disk',4,4); % structuring element r; lines 
closeBW = imclose(BWs,se); % close edges
figure (3), imshow(closeBW), title ('morphological closed edges');

BWdfill = imfill(closeBW, 'holes'); %fill holes
figure (4), imshow(BWdfill);
title('binary image with filled holes');

seD = strel('diamond',1); % structuring element
BWfinal = imerode(BWdfill,seD); % smoothen with imerode
%BWfinal = imerode(BWfinal,seD);
figure (5), imshow(BWfinal), title('smoothened image');

BWoutline = bwperim(BWdfill); % show outlines on original image
Segout = I; 
Segout(BWoutline) = 255; 
figure (6), imshow(Segout), title('outlines before segmentation on original image');

%% --- watershed transform based on intensities to seperate COD---

BWfinal_bw_unit8 = uint8 (BWdfill); %convert logical to integer
CaOxmask1 = I .* BWfinal_bw_unit8; 
figure (7), imshow (CaOxmask1), title ('masked image');

seD = strel('disk',5); % structuring element, decrease for higher segmentation
CaOxmask2 = imopen(CaOxmask1,seD); % smoothen with imopen
figure (8), imshow (CaOxmask2);
CaOxmask3 = imregionalmax (CaOxmask2);% create seeds for intensity based watershed
figure (9), imshow (CaOxmask3);


T = graydist (I, CaOxmask3); % distance transform of greyscale image, using seeds
figure (10), imshow (T); title ('distance transform of greyscale image');
Ld = watershed (T); % seeded watershed of distance transform
figure (11), imshow (Ld), title ('intensity based watershed');
figure (12), imshowpair (BWdfill, Ld, 'blend'), title ('intensity based watershed on binary image');
bwCOD = BWdfill; % add watershed segmentation to binary image
bwCOD (Ld == 0) = 0;
figure (13), imshow (bwCOD), title ('binary with intensity based segmentation');


%% --- watershed transform based on shape to seperate COM---

D = -bwdist(~BWdfill); % distance transform of binary image
mask = imextendedmin(D,1); % extended minima transform of h-minima transform (suppresses all minima <h)
figure (14), imshowpair(BWdfill,mask,'blend'), title ('minima imposition');
D2 = imimposemin(D,mask);% % modifies intensity image D; only regional minima where mask is not zero
figure (15), imshowpair (BWdfill, D2, 'blend');
Ld2 = watershed(D2); % shape based watershed segmentation to split touching COM
bwCOM = bwCOD; 
bwCOM(Ld2 == 0) = 0; %overlay shape based watershed on binary image
figure (16), imshow(bwCOM), title ('watershed seperated objects COM + COD');

BWoutline = bwperim(bwCOM);
Segout = I; 
Segout(BWoutline) = 255; 
figure (17), imshow(Segout), title('seperated object outlines on original image');

[B,L, N] = bwboundaries(bwCOM,'noholes');
figure (18) , imshow(label2rgb(L, @jet, [.5 .5 .5])), title ('seperated objects');
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end


%% --- extract features ---

bwBorder = imclearborder (bwCOM); % discard objects touching the border
stats = regionprops ('table', bwBorder, 'area','centroid', 'perimeter', 'boundingbox', 'eccentricity', 'majoraxislength', 'minoraxislength'); %extract shape measurments
centroid = cat (1,stats.Centroid); % reformat data to array
area = cat (1, stats.Area); % reformat the area data to array
perimeter = cat (1, stats.Perimeter); 
circularity = rdivide (area, perimeter);
majorAxis = cat(1, stats.MajorAxisLength);
minorAxis = cat (1, stats.MinorAxisLength);
eccentricity = cat (1, stats.Eccentricity);
boundingBox = cat (4, stats.BoundingBox);

statsGrayscale = regionprops ( bwBorder, I, {'MinIntensity' 'MaxIntensity' 'MeanIntensity' 'PixelValues'}); %extract intensity measurments
MinIntensity = cat (1, statsGrayscale.MinIntensity);
MinIntensity = double (MinIntensity);
MaxIntensity = cat (1, statsGrayscale.MaxIntensity);
MaxIntensity = double(MaxIntensity);
MeanIntensity = cat (1, statsGrayscale.MeanIntensity);

IntDist = MaxIntensity - MinIntensity; % calculate range of intensities
labels = bwlabel (bwBorder); 
count = numel (area);

s = regionprops (bwBorder, I, {'PixelValues'}); % calculate SD and variance of intensities
for k = 1:count
    s(k).StandardDeviation = std(double(s(k).PixelValues));
end
SD = cat (1, s.StandardDeviation);

for k= 1:count
    s(k).Variance = var(double(s(k).PixelValues));
end
 variance = cat (1, s.Variance);
   

%% ---extract texture features --- 
bwCOM_bw_unit8 = uint8 (bwBorder); %convert logical to integer
CaOxmaskFinal = I .* bwCOM_bw_unit8; 
figure (20), imshow (CaOxmaskFinal), title ('CaOx mask final');

for k= 1:count % crop image to subimage/particle and store in cell array
   singleObjects {k} = imcrop (CaOxmaskFinal, boundingBox (k,:));
end


for k= 1:count % calculate graycomatrix of each particle
   glcms1 {k} = graycomatrix (singleObjects{k});
end


for k = 1:count % exctract graycoprops
    contrast {k} = graycoprops(glcms1{k}, 'Contrast');
    correlation {k} = graycoprops(glcms1{k}, 'Correlation');
    energy {k} = graycoprops(glcms1{k}, 'Energy');
    homogeneity {k} = graycoprops(glcms1{k}, 'Homogeneity');
end


contrast1 = cell2mat(contrast); % reformat cell array to double
contrast2 = cat (1, contrast1.Contrast);
correlation1 = cell2mat(correlation); % reformat cell array to double
correlation2 = cat (1, correlation1.Correlation);
energy1 = cell2mat(energy); % reformat cell array to double
energy2 = cat (1, energy1.Energy);
homogeneity1 = cell2mat(homogeneity); % reformat cell array to double
homogeneity2 = cat (1, homogeneity1.Homogeneity);


%% --- feature descriptors to distinguish COM and COD ---

% choose one descriptor that distinguishes COM, COD and unclassified in training image
% e.g. MaxInt; variance, size - use to classify COM/COD in training images
% feature has to be adapted to each training image


CaOx = area >5; % discard noise --> all CaOx
CaOxVariance = variance(CaOx);
CaOxMeanInt = MeanIntensity(CaOx);
CaOxMinInt = MinIntensity(CaOx);
CaOxMaxInt = MaxIntensity(CaOx);
CaOxIntDist = IntDist (CaOx);
CaOxArea = area (CaOx);
CaOxPerimeter = perimeter (CaOx);
CaOxMajorAxis = majorAxis (CaOx);
CaOxMinorAxis= minorAxis (CaOx);
CaOxEccentricity = eccentricity (CaOx);
CaOxcircularity = circularity (CaOx);
CaOxcorrelation = correlation2 (CaOx);
CaOxcontrast = contrast2 (CaOx);
CaOxhomogeneity = homogeneity2 (CaOx);
CaOxenergy = energy2 (CaOx);

CaOxcentroid(:,1) = centroid (CaOx,1);
CaOxcentroid(:,2) = centroid (CaOx,2);



figure (21), imshow (Segout), title ('original image with MaxIntensity');
hold on;
 text(centroid(:,1), centroid (:,2), num2str(area, 4), 'color', 'red', 'fontsize',8); % add feature to distinguish COM/COD to outlined image to determine value
 
figure (22), imshow (Segout), title ('original image with MaxIntensity');
hold on;
 text(CaOxcentroid(:,1), CaOxcentroid (:,2), num2str(CaOxArea, 4), 'color', 'red', 'fontsize',8); % add feature to distinguish COM/COD to outlined image to determine value
 

class = CaOxArea <500; %classify COM from CaOx by choosing one feature, based on input image, e.g.variance


class1 = double (class);



%centroidCOM(:,1) = CaOxcentroid(class1); % new array with localization info (centroids, column 1)
%centroidCOM(:,2) = CaOxcentroid(class1,2); %new array with localization info (centroids, column 2)
%COMMaxInt = MaxIntensity (COM);
%COMarea = CaOxArea (class);
%varCOM = varianceariance (COM)

%figure (23), imshow (Segout), title ('original image with MaxIntensity');
%hold on;
 %text(centroidCOM(:,1), centroidCOM (:,2), num2str(class1, 4), 'color', 'red', 'fontsize',8); % add feature to distinguish COM/COD to outlined image to determine value
COM = class1 == 1;
COD = class1 == 0;
rods = class1 == 2;
dirt = class1 == 3;

[B1,L1,A1] = bwboundaries (bwBorder);
B2 = B1(CaOx == 1);
Bcom = B2(class1 == 1);
Bcod = B2 (class1 == 0);
Brod = B2 (class1 == 2);
Bdirt = B2 (class1 ==3);


figure (22) , imshow(I), title ('classified objects'); % plot classification of segmented crystals on original image
hold on
text(CaOxcentroid(:,1), CaOxcentroid (:,2), num2str(class1,1), 'color', 'black', 'fontsize',10);
hold on 
for k = 1:length(Bcom)
   boundary = Bcom{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

for k = 1:length(Bcod)
   boundary = Bcod{k};
   plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
end

for k = 1:length(Brod)
   boundary = Brod{k};
   plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2)
end

for k = 1:length(Bdirt)
   boundary = Bdirt{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% --- create categorial array with features ---
% arrange feature descriptors with COM/COD classification in table
% this table (combined from all training images) is used as the input for
% supervised classification (use classification learner app)
featuresAll = [class1 CaOxVariance, CaOxMeanInt, CaOxMinInt, CaOxMaxInt, CaOxIntDist, CaOxArea, CaOxPerimeter, CaOxEccentricity, CaOxMajorAxis, CaOxMinorAxis, CaOxcircularity, CaOxcorrelation, CaOxhomogeneity, CaOxenergy, CaOxcontrast]; % COM distinguishes if COM or COD; 
featuresAll1 = array2table (featuresAll, 'VariableNames', {'class' 'variance' 'MeanInt' 'MinInt' 'MaxInt' 'IntDis' 'area' 'perimeter' 'eccentricity' 'majorAxis' 'minorAxis' 'circularity' 'correlation' 'homogeneity' 'energy' 'contrast'});
writetable (featuresAll1, 'Classification09_exp09_NaOx_w205_train.csv');
