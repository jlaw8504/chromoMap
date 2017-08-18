pathname=uigetdir('C:\Users\lawrimor\Documents\MATLAB\ChromoMap\5000trimmed_MSD_analysis\','Pick the files');
fileList = dir(pathname);
for n = 1:64
    fileName = strcat(pathname,'\','coord_summ',num2str(n),'.csv');
    csvMatrix(:,:,n) = csvread(fileName,1,0);
end
[csvLength, ~, ~] = size(csvMatrix);
%% Check if all X's are all positive or all negative and take abs(x) if no point crosses midzone
%remember that X axis is the Z axis in the original chromoShake (column 4)
absTF = 1;
for j = 1:64
    for k = 1:csvLength
        signTF = sign(csvMatrix(k,4,j));
        if signTF ~= 0
            arrayTF(k) = signTF;
        else
            arrayTF(k) = 1;
            disp('I found a 0 value');
        end
    end
    sumArray = sum(arrayTF);
    sumArray = abs(sumArray);
    if sumArray ~= csvLength
        absTF = 0;
    end
end
if absTF == 1
   disp('No strand crosses the midpoint, can use abs(x) values.');
   csvMatrix(:,4,:) = abs(csvMatrix(:,4,:)); %abs of x for heatmap since model is symmetrical.
else
    disp('These points cross the midpoint, cannot use abs(x) values! Subtracting 635 nm Pole Position');

    for n = 1:64
        if mod(n,2) == 1 % if odd
            csvMatrix(:,4,n) = 6.35e-7 - csvMatrix(:,4,n); % odd numbered coord_summ X coords are positive at beginning of simulation
        else
            csvMatrix(:,4,n) = csvMatrix(:,4,n) + 6.35e-7; % even numbered coord_summ X coords are negative at beginning of simulation
        end
    end
end
%% Trim data to user specified offset
prompt = ['Your files have',' ', num2str(csvLength),' ', 'timepoints. How many timepoints should we disregard?\n'];
offset = input(prompt);
% In the simulation XYZ correlates to ZYX for what we could see under the
% microscope
x3D = csvMatrix((offset+1):end,4,:);
y3D = csvMatrix((offset+1):end,3,:);
z3D = csvMatrix((offset+1):end,2,:);
permX = permute(x3D,[1 3 2]);
X = reshape(permX,[],size(x3D,2),1);
permY = permute(y3D,[1 3 2]);
Y = reshape(permY,[],size(y3D,2),1);
permZ = permute(z3D,[1 3 2]);
Z = reshape(permZ,[],size(z3D,2),1);
%% Trim data to remove limit arrays in Z (the new Z we just permuted)
XYZ = horzcat(X,Y,Z);
[XYZ_length,~] = size(XYZ);
counter = 1;
Zmin = min(Z) * 10^9;
Zmax = max(Z) * 10^9;
disp(strcat('The range in Z (in nm) is',num2str(floor(Zmin)),':',num2str(floor(Zmax))));
disp('What range would you like to specify for Z?');
prompt_Zmin = ['Z minimum in nanometers:'];
min_bound = input(prompt_Zmin);
prompt_Zmax = ['Z maximum in nanometers:'];
max_bound = input(prompt_Zmax);
for n = 1:XYZ_length
    Z_check = XYZ(n,3);
    Y_check = XYZ(n,2);
    if Z_check >= (min_bound/10^9) && Z_check <= (max_bound/10^9)
        if Y_check >= (0/10^9) || Y_check <= (0/10^9)
            XYZ_trim(counter,:) = XYZ(n,:);
            counter = counter + 1;
        end
    end
end
X = XYZ_trim(:,1);
Y = XYZ_trim(:,2);
Z = XYZ_trim(:,3);
%% Calculate the fraction of points proximal to spindle axis
distal_counter = 0;
for n = 1:length(XYZ_trim(:,2))
if XYZ_trim(n,2) <= (-60*10^-9) || XYZ_trim(n,2) >= (60*10^-9)
distal_counter = distal_counter +1;
end
end

fraction_prox = (1 - distal_counter/length(XYZ_trim(:,2)))
%% Calculate Average Radial Distance
for m = 1:length(Z)
    radialDist(m) = sqrt(Z(m)^2 + Y(m)^2);
end
avgRadialDist_nm = (mean(radialDist)) * 10^9
stdRadialDist_nm = std(radialDist) * 10^9

%% Calculate Std in X after mirroring with abs(x), only possible for X coords that do not switch sign
if absTF == 1
    absX = abs(X);
    avgAbsX_nm = (mean(absX)) * 10^9
    stdAbsX_nm = (std(absX)) * 10^9
end

%% Create a Heatmap index
%setting up the X dimension on the array size by 64.8 nm pixels
%this is set to match the size of the images used for the 6.8 kb array
%in Stephens et al., 2011
Y = abs(Y); %Take the abs values in preparation for mirroring the H matrix
if absTF == 1
    X = 6.35e-7 - X; % put coordinates relative to a pole at 635 nm (AVG spindle length of Azide cells)
end
%Calculate mean and Std distance in X from pole
avgPoleX_nm = (mean(X)) * 10^9
stdPoleX_nm = (std(X)) * 10^9
if absTF == 1
    xDim = linspace((-1.296e-7),(7.776e-7),15)';
    yDim = linspace((-6.48e-7),(6.48e-7),21)';
else
    xDim = linspace((-1.296e-7),(1.3608e-6),24)';
    yDim = linspace((-6.48e-7),(6.48e-7),21)';
end
H = zeros(length(yDim),length(xDim));
for l = 1:length(X)
    countX = dsearchn(xDim, X(l));
    countY = dsearchn(yDim, Y(l));
    H(countY,countX) = H(countY,countX) + 1 ;
end
for n = 1:10
    H(n,:) = H((22-n),:);
end
%% Create the heatmap
% code taken from Matthew Larson's Heatmap make program 
data_name = 'H';
activedata=eval(data_name);
interpnumber = 2;
activeinterp=interp2(activedata,interpnumber);%linear interpolate data
activeinterp=activeinterp/max(max(activeinterp));%standardize to max=100%
pixelsize=input('What is your initial bin size in nanometers?\n'); %input pixel size
plottitle=input('What is the title of you plot?\n','s'); %input title
eval([data_name 'interp' num2str(interpnumber) '= activeinterp;'])%create name for interpolated data
eval(['interpvarname =''' data_name 'interp' num2str(interpnumber) ''';'])
figure
%create new plot
imagesc(activeinterp)
title(plottitle)
xlabel('Distance (nm)')
ylabel('Distance (nm)')
%adjust axis labels
xlabels=round([0:2:(size(activedata,2)-1)]*pixelsize);
xlabels=(mat2cell(xlabels,1,ones(1,length(xlabels))));
ylabels=round(fliplr(([0:2:(size(activedata,1)-1)]-floor(size(activedata,1)/2))*pixelsize));
ylabels=(mat2cell(ylabels,1,ones(1,length(ylabels))));
set(gca,'XTick',1:(2*2^interpnumber):size(activeinterp,2))
set(gca,'XTickLabel',xlabels)
set(gca,'YTick',1:(2*2^interpnumber):size(activeinterp,1))
set(gca,'YTickLabel',ylabels)
colorbar
colormap hot
%colormap jet %remove % at the beginning of this line for jet (rainbow) colormap
axis image
