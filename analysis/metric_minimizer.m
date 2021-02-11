
close all;

fileName           = 'test23';  % No file extension
chiSquareMapPlotOn = 1;
SNR = 25;

smoothingFactors = 1:1:2;
decodeVals       = -1:0.1:1;
decodeVals       = 1;

chiSquareArray = zeros(size(smoothingFactors, 2), size(decodeVals, 2));
counter1 = 1;

% Retrieve original drawn image PDF for comparison
[z, fig] = getOriginalImageData(fileName);

for smoothingFactor = smoothingFactors
    counter2 = 1;
    for decodeVal = decodeVals

        % Image reconstruction algorithm
        [image_estimate, RL_PSFest] = reconstructImageMethod1(fileName, SNR, [50 300], smoothingFactor, decodeVal);

        % Resize image to cast on to 256x256 image 
        image_estimate = resizem(image_estimate, 256/size(image_estimate,1));
        
        figure(3); subplot(1,2,1); cla; hold on;
        contourf(z); title("Image Reconstruction");
        contour(image_estimate); hold off;
        
        errorFigNum = 3;
        chiSquareArray(counter1, counter2) = computeChiSquared(image_estimate, z, chiSquareMapPlotOn, errorFigNum);
        
        counter2 = counter2 + 1;
    end
    counter1 = counter1 + 1;
end

[x,y] = meshgrid(smoothingFactors, decodeVals);

figure();
surf(x, y, chiSquareArray');
xlabel("Smoothing Factors");
ylabel("Decoder Null Values");
zlabel("\chi^2 Score");


function [image_est, PSF_est] = reconstructImageMethod1(fileName, SNR, energyRange, smoothingFactor, decodeVal)

lowE  = energyRange(1);
highE = energyRange(2);

darkDetectorID = 14;

hitFileName    = sprintf("../data/hit_%s.csv", fileName);    
signalFileName = sprintf("../data/signal_%s.csv", fileName);

data = importfile_resultsFile(signalFileName);
hits = importfile_resultsFile(hitFileName);

detectorIDs = unique([data.det; hits.det]);

Nsig       = size(data, 1); % c/s/intrument, total signal photons 
Nintrinsic = 20; % c/s/det, intrinsic detector noise
dt         = 10; % sec, integration time
ndet       = 11; % number of detectors

perDetRate = Nsig^2 * dt / SNR^2 - Nsig - Nintrinsic * ndet;

if perDetRate < 0; perDetRate = 0; end

perDetRate = 0;

pixelArray     = cell(length(detectorIDs), 2);
rawIm          = zeros(16,16);
darkIm         = zeros(16,16);
figure(1);
for detector = 1:length(detectorIDs)
    
    iPixels = [data.i(data.det == detectorIDs(detector) & data.E > lowE & data.E < highE) ; ...
        hits.i(hits.det == detectorIDs(detector) & hits.E > lowE & hits.E < highE)];
    
    jPixels = [data.j(data.det == detectorIDs(detector) & data.E > lowE & data.E < highE) ; ...
        hits.j(hits.det == detectorIDs(detector) & hits.E > lowE & hits.E < highE)];
   
    binnedData = histogram2(iPixels, jPixels, 0:1:16, 0:1:16);
    
    pixelArray{detector,1} = detectorIDs(detector);
    pixelArray{detector,2} = binnedData.Values;
    
    if detectorIDs(detector) ~= darkDetectorID 
        rawIm = rawIm + binnedData.Values + poissrnd(perDetRate/256, [16 16]);
    else
        darkIm = binnedData.Values + poissrnd(perDetRate/256, [16 16]);
    end
end
close 1;

% Normalize counts after image coaddition
%rawIm = rawIm / 11;

%load('CA_files/decoder.mat','decoder');
%load('CA_files/mask.mat'   ,'mask');
load('CA_files/NTHT_MURA_array_test.mat', 'mask','decoder');


% This value ensures that the sum over the image is equal to 0
% G = -1

% This value ensures that the sum over the image is equal to the counts
% over the detector. Lambda is the average height of the sidelobes of the
% autocorrelation of the mask with itself.
%
% G = -lambda / (M - lambda) = -60/(145-60)
%decoder(decoder == -1) = -0.7059;
%decoder(decoder == -1) = -0.1423; % NTHT

%decoder(decoder == -1) = decodeVal;

avgSigCounts = sum(sum(rawIm));
backgroundCounts  = sum(sum(darkIm));

totalFlux = 4.656 * (avgSigCounts - backgroundCounts);

naivePSF = conv2(mask, decoder);

%rawIm = rawIm - darkIm;

%rawIm = rawIm - mean(rawIm);
rawDeconv = conv2(rawIm, decoder);
%rawDeconv(rawDeconv < 0) = 0;

rawDeconv = interp2(rawDeconv, 1.1);
[image_est, PSF_est] = deconvblind(rawDeconv, naivePSF);

% Bypass iterative Lucy-Richardson deconvolution
%image_est = rawDeconv;
%PSF_est   = naivePSF;

%{
image_est = rawIm;
PSF_est   = naivePSF;
%}

% Normalize s.t. sum across image is 1
image_est = image_est / sum(sum(image_est));

% Assign signal flux to image and correct for inversion
image_est = rot90(image_est * totalFlux, 1);

image_est = imgaussfilt(image_est, smoothingFactor);

end

function [z, fig] = getOriginalImageData(fileName)

% Open original image and get data back from it
fig = openfig(['./originalImages/' fileName '.fig']);

dataObjs = findobj(fig,'-property','ZData');
z = resizem(dataObjs(2).ZData, 256/500);

close(fig);
end