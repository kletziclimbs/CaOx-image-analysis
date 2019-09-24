function results = myimfcn(I)
%function for batch processing app to segment CaOx 20x images and classify particles using the previously trainedClassifier09_4

%Image Processing Function
%
% I      - Input image.
% RESULTS - A scalar structure with the processing results.
%


%--------------------------------------------------------------------------
% Auto-generated by imageBatchProcessor App. 
%
% When used by the App, this function will be called for every input image
% file automatically. I contains the input image as a matrix. RESULTS is a
% scalar structure containing the results of this processing function.
%
%--------------------------------------------------------------------------


%input images - CaOx images 20x objective - COM/COD classification analysis 09 
%uses trained classifier 09_4

%output: count, mean and median size of all 4 classes: COD, COM, rods, dirt



%% ---detect edges---
Ismooth = imgaussfilt(I,1.4);
threshold = 0.03; % edge detection 
fudgeFactor = 0.5;
BWs = edge(Ismooth,'sobel', threshold * fudgeFactor); %make binary image of edges

%% --- morphologically close image and fill holes--
se = strel('disk',4,4); % structuring element r=2; lines =4
closeBW = imclose(BWs,se); % close edges

BWdfill = imfill(closeBW, 'holes'); %fill holes

%% --- watershed transform based on intensities to seperate COD---

BWfinal_bw_unit8 = uint8 (BWdfill); %convert logical to integer
CaOxmask1 = I .* BWfinal_bw_unit8; 

seD = strel('disk',5); % structuring element, decrease for higher segmentation
CaOxmask2 = imopen(CaOxmask1,seD); % smoothen with imopen
CaOxmask3 = imregionalmax (CaOxmask2);% create seeds for intensity based watershed

T = graydist (I, CaOxmask3); % distance transform of greyscale image, using seeds
Ld = watershed (T); % seeded watershed of distance transform
bwCOD = BWdfill; % add watershed segmentation to binary image
bwCOD (Ld == 0) = 0;



%% --- watershed transform based on shape to seperate COM---

D = -bwdist(~BWdfill); % distance transform of binary image
mask = imextendedmin(D,1); % extended minima transform of h-minima transform (suppresses all minima <h)
D2 = imimposemin(D,mask);% % modifies intensity image D; only regional minima where mask is not zero
Ld2 = watershed(D2); % shape based watershed segmentation to split touching COM
bwCOM = bwCOD; 
bwCOM(Ld2 == 0) = 0; %overlay shape based watershed on binary image

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


for k= 1:count
    s(k).Variance = var(double(s(k).PixelValues));
end
 variance = cat (1, s.Variance);
   

%% ---extract texture features --- 
bwCOM_bw_unit8 = uint8 (bwBorder); %convert logical to integer
CaOxmaskFinal = I .* bwCOM_bw_unit8; 

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


%% --- create categorial array with features ---

CaOx = area >10; % discard noise --> all CaOx
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



% arrange feature descriptors with COM/COD classification in table
% this table is used as the input for
% supervised classification using a trained classifier

featuresAll = [CaOxVariance, CaOxMeanInt, CaOxMinInt, CaOxMaxInt, CaOxIntDist, CaOxArea, CaOxPerimeter, CaOxEccentricity, CaOxMajorAxis, CaOxMinorAxis, CaOxcircularity, CaOxcorrelation, CaOxhomogeneity, CaOxenergy, CaOxcontrast]; % COM distinguishes if COM or COD; 
featuresAll1 = array2table (featuresAll, 'VariableNames', {'variance' 'MeanInt' 'MinInt' 'MaxInt' 'IntDis' 'area' 'perimeter' 'eccentricity' 'majorAxis' 'minorAxis' 'circularity' 'correlation' 'homogeneity' 'energy' 'contrast'});

%% -- classify ---
%classes are COD = 0, COM = 1; n.d./rods = 2; dirt = 3
load ('trainedClassifier09_4.mat');
yfit = trainedClassifier09_4.predictFcn(featuresAll1);

%% --- generate results ---
COM = yfit == 1;
COD = yfit == 0;
rods = yfit == 2;
dirt = yfit == 3;

COMarea = CaOxArea (COM);
CODarea = CaOxArea (COD);
rodarea = CaOxArea (rods);
dirtarea = CaOxArea (dirt);

medianCOM = median (COMarea);
meanCOM = mean (COMarea);
medianCOD = median (CODarea);
meanCOD = mean (CODarea);
countCOM = numel (COMarea);
countCOD = numel (CODarea);
countRod = numel (rodarea);
countDirt = numel (dirtarea);
medianRod = median (rodarea);
meanRod = mean (rodarea);
medianDirt = median (dirtarea);
meanDirt = mean (dirtarea);

results.medianCOM = medianCOM
results.medianCOD = medianCOD
results.meanCOM = meanCOM
results.meanCOD = meanCOD
results.meanRod = meanRod
results.medianRod = medianRod
results.meanDirt = meanDirt
results.medianDirt = medianDirt
results.countCOM = countCOM
results.countCOD = countCOD
results.countRod = countRod
results.countDirt = countDirt
%results.bwCOM = bwCOM







%--------------------------------------------------------------------------
