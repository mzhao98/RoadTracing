%% Fuzzy Logic Image Processing
% This example shows how to use Fuzzy Logic Toolbox(TM) software for image
% processing. Specifically, this example shows how to detect edges in an
% image.
%
% An edge is a boundary between two uniform regions. You can detect an edge
% by comparing the intensity of neighboring pixels. However, because
% uniform regions are not crisply defined, small intensity differences
% between two neighboring pixels do not always represent an edge. Instead,
% the intensity difference might represent a shading effect.
%
% The fuzzy logic approach for image processing allows you to use
% membership functions to define the degree to which a pixel belongs to an
% edge or a uniform region.
%
% Copyright 2013 The MathWorks, Inc. 

%% Import RGB Image and Convert to Grayscale
% Import the image into MATLAB.
Irgb = imread('terrain_pic_2.jpeg');

%%
% |Irgb| is a 384 x 512 x 3 |uint8| array. The three channels of |Irgb|
% (third array dimension) represent the red, green, and blue intensities of
% the image.
%%
% Convert |Irgb| to grayscale so that you can work with a 2-D array instead
% of a 3-D array. Use the standard NTSC conversion formula to calculate
% the effective luminance of each pixel.
Igray = 0.2989*Irgb(:,:,1)+0.5870*Irgb(:,:,2)+0.1140*Irgb(:,:,3);

figure
image(Igray,'CDataMapping','scaled');
colormap('gray')
title('Input Image in Grayscale')
%%
% Alternatively, you can use the |rgb2gray| function in the Image
% Processing Toolbox(TM) software to convert |Irgb| to grayscale.
%
%% Convert Image to Double-Precision Data
%
% The Fuzzy Logic Toolbox software operates on double-precision numbers
% only. So, convert |Igray|, a |uint8| array, to a |double| array.
I = double(Igray); 
%%
% Because |uint8| values are in the [0 2^8-1] range, all elements of |I|
% are in that range too. Scale |I| so that its elements are in the [0 1]
% range.
%
classType = class(Igray);                 
scalingFactor = double(intmax(classType));  
I = I/scalingFactor;      
%%
% Alternatively, you can use the |im2double| function in the Image
% Processing Toolbox software to convert |Igray| to a scaled,
% double-precision image.
%% Obtain Image Gradient
%
% The fuzzy logic edge-detection algorithm for this example relies on the
% image gradient to locate breaks in uniform regions. Calculate the image
% gradient along the _x_-axis and _y_-axis.
Gx = [-1 1];
Gy = Gx';
Ix = conv2(I,Gx,'same');
Iy = conv2(I,Gy,'same');

figure
image(Ix,'CDataMapping','scaled')
colormap('gray')
title('Ix')

figure
image(Iy,'CDataMapping','scaled')
colormap('gray')
title('Iy')
%% 
% |Gx| and |Gy| are simple gradient filters. You convolve |I| with |Gx|,
% using the |conv2| function, to obtain a matrix containing the _x_-axis
% gradients of |I|. The gradient values are in the [-1 1] range.
% Similarly, you convolve |I| with |Gy| to obtain the _y_-axis gradients of
% |I|. You can use other filters to obtain the image gradients, such as the
% Sobel operator or the Prewitt operator. For information about how you can
% filter an image using convolution, see
% <http://www.mathworks.com/help/images/designing-and-implementing-linear-filters-in-the-spatial-domain.html#f16-20755
% Convolution>.
%
% Alternatively, if you have the Image Processing Toolbox software, you can
% use the |imfilter|, |imgradientxy|, or |imgradient| functions to obtain
% the image gradients.
%% Define Fuzzy Inference System (FIS) for Edge Detection
% Create a fuzzy inference system (FIS) for edge detection, |edgeFIS|.
%
edgeFIS = newfis('edgeDetection');
%%
% Specify the image gradients, |Ix| and |Iy|, as the inputs of |edgeFIS|.
edgeFIS = addvar(edgeFIS,'input','Ix',[-1 1]);
edgeFIS = addvar(edgeFIS,'input','Iy',[-1 1]);
%%
% Specify a zero-mean Gaussian membership function for each input. If the
% gradient value for a pixel is |0|, then it belongs to the zero membership
% function with a degree of |1|.
sx = 0.1;
sy = 0.1;
edgeFIS = addmf(edgeFIS,'input',1,'zero','gaussmf',[sx 0]);
edgeFIS = addmf(edgeFIS,'input',2,'zero','gaussmf',[sy 0]);
%%
% |sx| and |sy| specify the standard deviation for the zero membership
% function for the |Ix| and |Iy| inputs. You can change the values of |sx|
% and |sy| to adjust the edge detector performance. Increasing the values
% makes the algorithm less sensitive to the edges in the image and
% decreases the intensity of the detected edges.
%%
% Specify the intensity of the edge-detected image as an output of
% |edgeFIS|.
edgeFIS = addvar(edgeFIS,'output','Iout',[0 1]); 
%%
% Specify the triangular membership functions, white and black, for |Iout|.
wa = 0.1;
wb = 1;
wc = 1;
ba = 0;
bb = 0;
bc = 0.7;
edgeFIS = addmf(edgeFIS,'output',1,'white','trimf',[wa wb wc]);
edgeFIS = addmf(edgeFIS,'output',1,'black','trimf',[ba bb bc]);
%%
% As you can with |sx| and |sy|, you can change the values of |wa|, |wb|,
% |wc|, |ba|, |bb|, and |bc| to adjust the edge detector performance. The
% triplets specify the start, peak, and end of the triangles of the
% membership functions. These parameters influence the intensity of the
% detected edges.
% 
%%
% Plot the membership functions of the inputs/outputs of |edgeFIS|.
figure
subplot(2,2,1)
plotmf(edgeFIS,'input',1)
title('Ix')
subplot(2,2,2)
plotmf(edgeFIS,'input',2)
title('Iy')
subplot(2,2,[3 4])
plotmf(edgeFIS,'output',1)
title('Iout')
%% Specify FIS Rules
% Add rules to make a pixel white if it belongs to a uniform region.
% Otherwise, make the pixel black.
r1 = 'If Ix is zero and Iy is zero then Iout is white';
r2 = 'If Ix is not zero or Iy is not zero then Iout is black';
r = char(r1,r2);
edgeFIS = parsrule(edgeFIS,r);
showrule(edgeFIS)
%% Evaluate FIS
% Evaluate the output of the edge detector for each row of pixels in |I|
% using corresponding rows of |Ix| and |Iy| as inputs.
Ieval = zeros(size(I));
for ii = 1:size(I,1)
    Ieval(ii,:) = evalfis([(Ix(ii,:));(Iy(ii,:));]',edgeFIS);
end
%% Plot Results
figure
image(I,'CDataMapping','scaled')
colormap('gray')
title('Original Grayscale Image')

figure
image(Ieval,'CDataMapping','scaled')
colormap('gray')
title('Edge Detection Using Fuzzy Logic')
%% Summary
% You detected the edges in an image using a FIS, comparing the gradient of
% every pixel in the _x_ and _y_ directions. If the gradient for a pixel is
% not zero, then the pixel belongs to an edge (black). You defined the
% gradient as zero using Gaussian membership functions for your FIS inputs.
%
