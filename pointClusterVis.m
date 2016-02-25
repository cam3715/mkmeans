function pointClusterVis(clustering,trialnum)
%POINTCLUSTERVIS plots polar and cylinderical visualizations of clustering
% 
% Input:
%     clustering             - struct
%             .k             - number of means
%             .idx           - cluster assignments (indicies) for each trial
%             .distances     - distance from each point to each center
%     trialnum               - trial to visualize (default 1)
%
% Author: Camden Glenn Bock
% 598 Bates College, Lewistion, ME 04240
% cbock@bates.edu, camdenbock@gmail.com
% http://www.camdenbock.com
% December 2015; Last Revision: 12/30/2015
%
%% Copyright (C) 2016  Camden Bock - GPL v. 3.0
%
% 'This program is liscensed under GPL v3.0'
% 'This program is modified from Ahmad Alsahaf`s package: ...
% amjams/mixedkmeans'.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% 'This program comes with ABSOLUTELY NO WARRANTY;' ...
% 'for details type view source. This is free software, and' ...
% 'you are welcome to redistribute it display under certain' ...
% 'conditions; see <http://www.gnu.org/licenses/>', ...
% ('Copyright (C) 2016  Camden Bock, Bates College'))



%------------- BEGIN CODE --------------

if nargin < 2
    trialnum = 1;
elseif nargin < 1
    display('not enough arguments')
end

idx = clustering.idx(:,trialnum);
distances = clustering.distance(:,:,trialnum);
rownum = 1:length(idx);

[idx,I] = sort(idx);
distances = distances(I,:);

k = numel(unique(idx));
X = zeros(length(rownum),k+1);
Y = zeros(length(rownum),k+1);
Z = zeros(length(rownum),k+1);
C = zeros(length(rownum),k+1);

%plot regular polygon
deltaAngle = 2*pi/k;
theta = 0:deltaAngle:2*pi;

radiusReg = zeros(1,k+1);
for i=1:k+1
    radiusReg(i) = max(max(distances));
end
figure
subplot(1,2,1)
if gpuDeviceCount > 0
    thetaG = gpuArray(theta);
    radiusRegG = gpuArray(radiusReg);
    polar(thetaG,radiusRegG,'b');
else
    polar(theta,radiusReg,'-b');
end
hold on
for i=1:length(rownum)
    distancePoint = distances(rownum(i),:);
    %plot datapoint point polygon
    radiusPoint = zeros(1,k+1);
    parfor j=1:k
        radiusPoint(j) = distancePoint(j);
    end
    radiusPoint(k+1) = radiusPoint(1);
    if gpuDeviceCount > 0
        thetaG = gpuArray(theta);
        radiusPointG = gpuArray(radiusPoint);
        polar(thetaG,radiusPointG,'-r');
    else
        polar(theta,radiusPoint,'-r');
    end
    [x,y] = pol2cart(theta, radiusPoint);
    X(i,:) = x;
    Y(i,:) = y;
    Z(i,:) = rownum(i);
    for j=1:k
        C(i,j) = idx(i);
    end
end
hold off
subplot(1,2,2)
if gpuDeviceCount > 0
    Xg = gpuArray(X);
    Yg = gpuArray(Y);
    Zg = gpuArray(Z);
    Cg = gpuArray(C);
    surf(Xg,Yg,Zg,Cg);
else
    surf(X,Y,Z,C);
end
end