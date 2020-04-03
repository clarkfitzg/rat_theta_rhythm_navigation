function [ThetaIndex] = Get_ThetaIndex(EntityID,clusterID,Binwidth,maxXaxis)
% this function gets the mean of the autocorrelation bin and plots the
% location of the mean on the graph

%chahe@tcd.ie
flag = 0;
if nargin < 2
    error('\nplease enter the EntityID and clusterID');
elseif nargin == 2
    [counts,binEdges] = PlotCorrelationHistograms(EntityID,clusterID);
elseif nargin == 3
    [counts,binEdges] = PlotCorrelationHistograms(EntityID,clusterID,Binwidth);
elseif nargin == 4
    [counts,binEdges] = PlotCorrelationHistograms(EntityID,clusterID,Binwidth,maxXaxis);
end
hold on
counts(isnan(counts)) = 0;
if max(counts) == 0 || sum(counts) <= 20
    ThetaIndex = 0;
    h = findobj('Name','Correlogram');
set(h,'visible','on');
set(h,'Name','ThetaIndex');
    return;
end
mid = find(binEdges == 0,1);

m = max(counts(mid-50/Binwidth:mid+50/Binwidth));
Indi = find(counts >= .999*m);
if max(Indi) >= mid+50/Binwidth
    Indi = Indi(find(Indi <= mid+50/Binwidth));
end

counts(min(Indi):max(Indi)) = m;
if m ~= 0
counts = counts./m;
end
options = optimset('TolFun',1e-10,'MaxIter',500);
% mean_of_spikes = sum(counts(mid:end))/2; % median not the mean
% mean_ind = find(counts(mid:end) <= mean_of_spikes);
% for I=1:length(counts)-mid
%     if sum(counts(mid:mid+I)) >= mean_of_spikes
%         mean_time = binEdges(mid+I);
%         break;
%     end
% end
% plot([mean_time mean_time], get(gca, 'YLim'), 'k:')
B = fir1(60,[.04*Binwidth]);
smoothed =  filtfilt(B,[1],counts);
%plot(binEdges,smoothed,'b');
%counts = fliplr(counts);
h = findobj('Name','Correlogram');
set(h,'visible','on');
set(h,'Name','ThetaIndex');
hold on;

EndInd = find(binEdges >= 700,1);
startInd = find(binEdges <= -700,1,'Last');
counts = counts(startInd:EndInd);

binEdges = binEdges(startInd:EndInd);
%counts(find(counts <= 20 )) = 0;
smoothed = smoothed(startInd:EndInd);
mid = find(binEdges == 0,1);
% estimate the first parameter which is b and tou1
%[x1,resnorm] = lsqcurvefit(@fun1,[counts(mid) 10],binEdges,counts);
% now estimate the second set of parameters  c and tou2 and reevalate the fist two
%[x2,resnorm] = lsqcurvefit(@fun2,[x1 8000 .01 1000],binEdges,counts);% note 10 10 is just the initial guess
% now the third set of parameters a and w = 2*pi*f
%[x3,resnorm] = lsqcurvefit(@fun3,[x2 50],binEdges,counts); % note the first guess .01 is very important
% now x3 is the estimated parameters with this order 
% x(1) = b % x(2) = tou1 % x(3) = a % x(4) = f % x(5) = c % x(6) = tou2  
%ThetaIndex = x3(3)/x3(1);
%F3 = fun3(x3,binEdges);
%plot(binEdges,F3,'r','lineWidth',2);
smoothed = smoothed(mid:end);
slope = diff(smoothed);

peak_ = zeros(size(smoothed));
for I =2:length(slope)
    if slope(I) > 0 && slope(I-1) <=0
        peak_(I+1) = -1;
    elseif slope(I) < 0 && slope(I-1) >=0
        peak_(I+1) = 1;
    end
end
binEdgesN  = binEdges(mid:end);
countsN = counts(mid:end);
pos = find(peak_>0);
neg = find(peak_ < 0);
if (~isempty(pos) & ~isempty(neg))
if binEdgesN(pos(1)) < 50
    pos = pos(2:end);
end
neg = neg(2:end);
end
%plot(binEdges,smoothed);
%plot(binEdges(pos),smoothed(pos),'g*');
%plot(binEdges(neg),smoothed(neg),'b*');
% now calculate the theta index
%smoothed = smoothed./smoothed(pos(1));
Lpos = length(pos);
Lneg = length(neg);
%figure()
if (~isempty(pos) & ~isempty(neg)) & Lneg > 3
    %ampdf = mean(counts(pos(1):end))-mean(smoothed(pos(1):end));
    %if abs(ampdf) >= 0.1
    %    smoothed = smoothed + ampdf;
    %end
    %counts = counts./smoothed(pos(1));
    %smoothed = smoothed./smoothed(pos(1));
    
    mean(countsN(pos(1):end));
   % plot(binEdges,counts,'r')
    %hold on
    %plot(binEdges,smoothed,'g')
    %plot([binEdges(1):max(binEdges)],mean(smoothed),'k:')
   % plot([binEdges(1):max(binEdges)],mean(counts),'b:')
    if Lpos ~= Lneg
        if Lpos > Lneg
            pos = pos(1:length(pos)-(Lpos-Lneg));
        elseif Lpos <Lneg
            neg = neg(1:length(pos));
        end
    end
    Lpos = length(pos);
Lneg = length(neg);
    if Lneg >= 3
    for I=2:Lneg
        decay_constant(I-1) = ((binEdgesN(neg(I-1)-1)-binEdgesN(neg(I)-1))/(log(smoothed(neg(I)-1)./smoothed(neg(I-1)-1))));
        %b = counts(neg(I)
    end
    else 
        decay_constant = 100;
    end
     decay_constant = real(decay_constant);
    av_decay_constant3 = mean(decay_constant);
    b = smoothed(neg-1)./(exp(-binEdgesN(neg-1)./av_decay_constant3));
    av_b = mean(b);
    if isnan(av_b)
        av_b = 0;
    end
    if av_decay_constant3 <= 0
        av_decay_constant3 = 250;
    end
    [x4,resnorm] = lsqcurvefit(@fun4,[av_b av_decay_constant3],binEdges,counts);
    %a = (smoothed(pos)./(2*exp(-binEdges(pos)./av_decay_constant)))-(b./2);
    %b = (smoothed(neg)./exp(-binEdgesN(neg)./av_decay_constant)+ smoothed(pos)./exp(-binEdgesN(pos)./av_decay_constant))./2;
    
    if Lpos >= 3
     for I=2:Lpos
        decay_constant2(I-1) = (binEdgesN(pos(I)-1)-binEdgesN(pos(I-1)-1))/(log((smoothed(pos(I-1)-1)-(av_b.*exp(-binEdgesN(pos(I-1)-1)./x4(2))))./(smoothed(pos(I)-1)-(av_b.*exp(-binEdgesN(pos(I)-1)./x4(2))))));
        %b = counts(neg(I)
        if I<= 3
        w(I-1) = 2*pi./(binEdgesN(pos(I)-1)-binEdgesN(pos(I-1)-1));
        end
    end
    else 
        decay_constant2 = 100;
    end
    decay_constant2 = real(decay_constant2);
    %a =  smoothed(pos)./exp(-binEdgesN(pos)./av_decay_constant) - av_b;
    av_decay_constant2 = mean(decay_constant2);
    a =  (smoothed(pos-1) - ((av_b*exp(-binEdgesN(pos-1)./x4(2)))))./(2*exp(-binEdgesN(pos-1)./av_decay_constant2));
    av_a = mean(a);
    %[x2,resnorm] = lsqcurvefit(@fun2,[av_b x1(2) av_a mean(w) x1(2)],binEdges,counts);% note 10 10 is just the initial guess

    %av_b = mean(b);%+mean(a);
    %w = [1:length(pos)].*2*pi./binEdgesN(pos-1);
    
    
    
    
    
    % a = smoothed(pos)-smoothed(neg);
    % b = smoothed(neg);
    ThetaIndex = mean(av_a/av_b);
    
else
    ThetaIndex  = 0
    av_a = 0;
    av_b = countsN(80);
    w = 0.01;
    av_decay_constant2 =250;
    flag = 1;
    av_decay_constant3 = 250;
    
end
countsN(isnan(countsN)) = 0;
counts(isnan(counts)) = 0;

av_b = real(av_b);av_decay_constant2 = real(av_decay_constant2);av_a  = real(av_a);
%options = OPTIMSET('MaxFunEvals', 100000,'MaxIter',10000);
%[x1,resnorm] = lsqcurvefit(@fun1,[x2(1) x2(2)],binEdges,counts);
% now estimate the second set of parameters  c and tou2 and reevalate the fist two
if min(w) == max(w)
    w(length(w)+1) = max(w)+max(w)*0.15;
    w(length(w)+1) = max(w)-max(w)*0.15;
end
[x1,resnorm] = lsqcurvefit(@fun1,[av_a mean(w) 10],binEdges,counts,[0 min(w) 0],[1 max(w) inf])
[x2,resnorm] = lsqcurvefit(@fun2,[av_b av_decay_constant3 x1],binEdges,counts,[0  -inf 0 min(w)],[1 inf  1 max(w)]);% note 10 10 is just the initial guess
% now the third set of parameters a and w = 2*pi*f
if counts(mid) ~= 0
%[x3,resnorm] = lsqcurvefit(@fun3,[counts(mid)-x2(1)-2*x2(3) 10 x2],binEdges,counts,[0 1],[counts(mid) 50]); % note the first guess .01 is very important
[x3,resnorm,resi,EXITFLAG] = lsqcurvefit(@fun3,[counts(mid)-x2(1)-2*x2(3) 10  x2],binEdges,counts,[ 0 0 0  -inf 0 min(w)],[counts(mid) 50 1 inf  1 max(w)],options);
% now x3 is the estimated parameters with this order 
% x(1) = b % x(2) = tou1 % x(3) = a % x(4) = f % x(5) = c % x(6) = tou2  

ThetaIndex = (x3(5)/x3(3))*100;
if x3(3) <=1e-5
    ThetaIndex = 0;end
ThetaIndex = round(ThetaIndex)/100;
F3 = fun3(x3,binEdges);
F3 = F3.*m;
plot(binEdges,F3,'r','lineWidth',2)

xiiii = 1;
% figure()
% plot(counts)
% hold on
% plot(F3./m,'r')
else 
    ThetaIndex = 0;
end
if ThetaIndex > 0.1 && flag ==1;
     ThetaIndex = 0;
end
%x3 = [av_a;mean(w);av_decay_constant;av_b];


    
