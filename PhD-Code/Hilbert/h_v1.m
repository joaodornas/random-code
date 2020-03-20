function h_v1

dataLength = 900000;

% %% Load DataSet Sanity
% data = load('sanity.mat');
% X = data.Di(1:dataLength);
% Y = data.Ci(1:dataLength);
% Trans = data.Tevent(data.Tevent<dataLength);
% dataLabel = 'Sanity';

%% Load DataSet Robin
data = load('RobinTimeSeries.mat');
trans = load('Time-Transitions.mat');
X = data.X2dt(1:dataLength) - data.Y2dt(1:dataLength); %% Decision
Y = data.X1dt(1:dataLength) - data.Y1dt(1:dataLength); %% Evidence
Trans = trans.trans.time(trans.trans.time<dataLength);
dataLabel = 'Robins Model';

% figure
% plot(1:dataLength,X,'b');
% hold on
% plot(1:dataLength,Y,'r');

%% Hilbert Transform
[xr, xw, xwt, hX] = hilbertTransform(X);
[yr, yw, ywt, hY] = hilbertTransform(Y);

xw_trans = xw(Trans);
yw_trans = yw(Trans);

xwt_positive = xwt;
ywt_positive = ywt;

xwt_positive(xwt_positive<0) = 0;
ywt_positive(ywt_positive<0) = 0;

% xw_seg = getSegmentsTrans(xw,Trans,nTrans,segmentsLength);
% yw_seg = getSegmentsTrans(yw,Trans,nTrans,segmentsLength);

% %% Granger Causality
% alpha = 0.05;
% max_lag = 5;
% 
% [xw_yw_F, xw_yw_CV] = granger_cause(yw,xw,alpha,max_lag);
% [yw_xw_F, yw_xw_CV] = granger_cause(xw,yw,alpha,max_lag);
% 
% 
% %% Plots
% 
% % Granger Causality between Hilbert Transform Phases with Transitions
% disp(strcat('Decision Granger-cause Evidence:',int2str(xw_yw_F>xw_yw_CV)));
% disp(strcat('F-value:',num2str(xw_yw_F)));
% disp(strcat('Critical Value:',num2str(xw_yw_CV)));
% 
% disp(strcat('Evidence Granger-cause Decision:',int2str(yw_xw_F>yw_xw_CV)));
% disp(strcat('F-value:',num2str(yw_xw_F)));
% disp(strcat('Critical Value:',num2str(yw_xw_CV)));

% Hilbert Transform Orbits with Transitions

% figure
% plot(xw,yw,'b');
% hold on
% plot(abs(xw_trans),abs(yw_trans),'ro');
% title(strcat('Orbits (',dataLabel,')'));
% xlabel('Decision Phases');
% ylabel('Evidence Phases');
% 
%Hilbert Transform Amplitude Phase Frequency Transitions

vector_size = 500000;

limit = 0;

while (limit < 1) || (limit < vector_size)
    
    limit = ceil(rand(1)*1000)*dataLength/1000;
    
end

xr = xr((limit-vector_size+1):limit);
yr = yr((limit-vector_size+1):limit);
xw = xw((limit-vector_size+1):limit);
yw = yw((limit-vector_size+1):limit);
xwt = xwt((limit-vector_size):limit-1);
ywt = ywt((limit-vector_size):limit-1);
xwt_positive = xwt_positive((limit-vector_size):limit-1);
ywt_positive = ywt_positive((limit-vector_size):limit-1);

sel_trans = intersect(Trans,(limit-vector_size+1):limit);
sel_trans = sel_trans - (limit-vector_size+1);

% high_freq = find(abs(ywt(1:end))>1);
% 
% yw(high_freq) = yw(high_freq+1);

% xw = smooth(xw);
% yw = smooth(yw,100);
      
% figure
% plot(1:vector_size,xr,'b');
% hold on
% plot(1:vector_size,xw,'r');
% hold on
% plot(1:length(xwt_positive),xwt_positive,'g');
% hold on
% plot(sel_trans,ones(length(sel_trans),1),'ko');
% legend({'Amplitude' 'Phase' 'Frequency' 'Transitions'});
% title(strcat('Decision (',dataLabel,')'));
% 
% figure
% plot(1:vector_size,yr,'b');
% hold on
% plot(1:vector_size,yw,'r');
% hold on
% plot(1:length(ywt_positive),ywt_positive,'g');
% hold on
% plot(sel_trans,ones(length(sel_trans),1),'ko');
% legend({'Amplitude' 'Phase' 'Frequency' 'Transitions'});
% title(strcat('Evidence (',dataLabel,')'));
% 
% figure
% plot(1:vector_size,yw,'b');
% hold on
% plot(1:vector_size,xw,'r');
% hold on
% plot(sel_trans,ones(length(sel_trans),1),'ko');
% legend({'Evidence' 'Decision''Transitions'});
% title(strcat('Phases (',dataLabel,')'))

% 
% figure
% plot(1:dataLength,xw,'b');
% hold on
% plot(1:dataLength,yw,'r');
% hold on
% plot(Trans,ones(length(Trans),1),'ko');

% 
% figure
% plot(1:dataLength,imag(hX),'b');
% hold on
% plot(1:dataLength,imag(hY),'r');
% hold on
% plot(Trans,ones(nTrans,1),'ko');
% legend({'Decision Hilbert' 'Evidence Hilbert' 'Transitions'});
% 


%% Hilbert Transform Frequencies
% % yr_norm = yr./max(yr);
% % xr_norm = xr./max(xr);
% 
% figure
% % plot(1:dataLength,yr_norm,'b--');
% % hold on
% % plot(1:dataLength,xr_norm,'r--');
% % hold on
% plot(1:length(xwt_positive),xwt_positive,'r');
% hold on
% plot(1:length(ywt_positive),ywt_positive,'b');
% hold on
% plot(Trans,ones(length(Trans),1)-0.5,'ko');
% % legend({'Evidence Aplitude (norm.)' 'Decision Aplitude (norm.)' 'Evidence Frequency' 'Decision Frequency'});
% legend({'Evidence Frequency' 'Decision Frequency' 'Transitions'});
% title(dataLabel);

% %% Phases per Segment
% 
% offset = -segmentsLength:segmentsLength;
% 
% figure
% for k=1:size(xw_seg,1)
%     
%     plot(offset,abs(xw_seg(k,:)),'b');
%     hold on
%     
% end
% legend({'Decision Phases Around Transitions'});
% 
% figure
% for k=1:size(yw_seg,1)
%     
%     plot(offset,abs(yw_seg(k,:)),'r');
%     hold on
%     
% end
% legend({'Evidence Phases Around Transitions'});
% 
% figure
% for k=1:size(xw_seg,1)
%     
%     plot(abs(xw_seg(k,:)),abs(yw_seg(k,:)),'k');
%     hold on
%     
% end
% title({'Phases Around Transitions'});
% xlabel('Decision');
% ylabel('Evidence');

% [new_xw, new_yw] = takeOutJumps(xw,yw,Trans);

% color{1} = 'b';
% color{2} = 'k';
% color{3} = 'g';
% color{4} = 'm';
% color{5} = 'y';
% color{6} = 'c';
% color{7} = 'r';
% 
% for h=1:7
%   
%     phase_x = h*10;
%     phase_y = 0;
% 
%     nTrans = 10;
% 
%     [H(h).seg_x, H(h).seg_y, H(h).trans_x_angle_idx, H(h).trans_y_angle_idx, H(h).phase_trans_x, H(h).phase_trans_y, H(h).phase_offset] = getSegmentsPhases(xw,yw,xwt,ywt,sel_trans,nTrans,phase_x,phase_y);
% 
% end
% 
% for h=1:7
%     
%     for g=1:size(H(h).seg_x,2)
% 
%         plot(H(h).seg_y(g).phases(H(h).trans_x_angle_idx(g)),H(h).seg_x(g).phases(H(h).trans_x_angle_idx(g)),'r.','MarkerSize',30);
%         hold on
%         
%         plot(1:length(H(h).seg_x(g).phases),H(h).seg_x(g).phases,color{g});
% 
%            min_y = min(H(h).seg_y(g).phases);
%            min_x = min(H(h).seg_x(g).phases);
%     
%            max_y = max(H(h).seg_y(g).phases);
%            max_x = max(H(h).seg_x(g).phases);
%     
%            min_h_phases(h,g) = min([min_y min_x]);
%            max_h_phases(h,g) = max([max_y max_x]);
%     
%     end
%     
%     min_phases = min(min_h_phases);
%     max_phases = max(max_h_phases);
%     
% end
% 
% title('Orbit');
% ylim([min(min_phases)-10 max(max_phases)+10]);
% xlabel('Phases - Evidence');
% ylabel('Phases - Decision');

phase_x = 0;
phase_y = 90;

nTrans = 5;

[seg_x, seg_y, trans_x_angle_idx, trans_y_angle_idx, phase_trans_x, phase_trans_y, phase_offset] = getSegmentsPhases(xw,yw,xwt,ywt,sel_trans,nTrans,phase_x,phase_y);

color{1} = 'b';
color{2} = 'k';
color{3} = 'g';
color{4} = 'm';
color{5} = 'y';

% for g=1:size(seg_x,2)
%     
%     l_y(g) = length(seg_y(g).phases);
%     l_x(g) = length(seg_x(g).phases);
%     
% end
% 
% min_l_y = min(l_y);
% min_l_x = min(l_x);
% 
% min_l = min([min_l_y min_l_x]);
% 
% for g=1:size(seg_x,2)
%     
%    seg_x(g).phases = seg_x(g).phases(1:min_l);
%    seg_y(g).phases = seg_y(g).phases(1:min_l);
%     
% end

timePoints = 5;

for g=1:size(seg_x,2)
   
    start = trans_x_angle_idx(g) - timePoints - 1;
    end_ = trans_x_angle_idx(g) + timePoints;
    
    seg_x(g).phases = seg_x(g).phases(start:end_-1);
    seg_y(g).phases = seg_y(g).phases(start:end_-1);
    
    trans_x_angle_idx(g) = timePoints + 1;
    trans_y_angle_idx(g) = timePoints + 1;
    
end

figure
for g=1:size(seg_x,2)
     
    phases(g).x = seg_x(g).phases - seg_x(g).phases(trans_x_angle_idx(g));
    phases(g).y = seg_x(g).phases - seg_y(g).phases;
    
    plot(phases(g).y,phases(g).x,color{g});
    
    hold on
     
    min_y = min(seg_y(g).phases);
    min_x = min(seg_x(g).phases);
    
    max_y = max(seg_y(g).phases);
    max_x = max(seg_x(g).phases);
    
    min_phases(g) = min([min_y min_x]) - 1;
    max_phases(g) = max([max_y max_x]) + 1;
    
end

for g=1:size(seg_x,2)

    plot(phases(g).y(trans_y_angle_idx(g)),phases(g).x(trans_x_angle_idx(g)),'r.','MarkerSize',30);
    hold on
    
end
title('Orbit');
%ylim([min(min_phases) max(max_phases)]);
xlabel('Phases - Evidence');
ylabel('Phases - Decision');

figure
for g=1:size(seg_x,2)
    
    plot(1:length(seg_y(g).phases),seg_y(g).phases,color{g});
    hold on
    
    plot(trans_y_angle_idx(g),seg_y(g).phases(trans_y_angle_idx(g)),'r.','MarkerSize',30);
    hold on
    
end
title('Evidence');
ylim([min(min_phases) max(max_phases)]);
xlabel('Time Points');
ylabel('Phases');


figure
for g=1:size(seg_x,2)
    
    plot(1:length(seg_x(g).phases),seg_x(g).phases,color{g});
    hold on
    
    plot(trans_x_angle_idx(g),seg_x(g).phases(trans_x_angle_idx(g)),'r.','MarkerSize',30);
    hold on
    
end
title('Decision');
ylim([min(min_phases) max(max_phases)]);
xlabel('Time Points');
ylabel('Phases');


% xlim([min(seg_x) max(seg_x)]);
% ylim([min(seg_y) max(seg_y)]);

 function [amp, phase, freq, hW] = hilbertTransform(W)
        
       hW = hilbert(W);
     
       amp = abs(hW);
       phase = angle(hW);
       freq = diff(phase); 
        
 end

% function [uni_wx, uni_wy, dist] = condPhaseDis(wx,wy,tau)
%     
%     uni_wy = sort(unique(wy));
%     
%     uni_wx = sort(unique(wx));
%     
%     dist = zeros(length(uni_wy),length(uni_wx));
%     
%     for i=1:length(uni_wy)
%         
%         idx = find(wy==uni_wy(i));
%         
%         post_wx = wx(idx+tau);
%         
%         [uni_post_wx, numUni_post_wx] = count_unique(post_wx);
%         
%         numUni_post_wx = numUni_post_wx./max(numUni_post_wx);
%         
%         
%         
%         dist(i,:
%         
%     end
%         
%  end

function [new_x, new_y] = takeOutJumps(x,y,Transitions)
        
    t = 1;
    new = 0;
    last = 0;
    for t=1:(length(x)-1)
        
        new = x(t);

        i = 0;
        
        if (new<0) && (last>0)
            
            i = i + 1;
            
            inSeg = 1;
            k = 1;
            j = 1;
            while inSeg == 1
                
                new_y(i).seg(k) = y(t);
                new_x(i).seg(k) = x(t);
                
                if find(Transitions==t) == 1
                    
                    new_y(i).trans(j) = k;
                    new_x(i).trans(j) = k;
                    
                    j = j + 1;
                    
                end
                
                k = k + 1;
                
                if (x(t) > 0) && (x(t+1) < 0)
                    
                    inSeg = 0;
                    
                end
                
            end
            
        end
        
        last = new;
        
    end
            
            
end     
            
function [seg_x, seg_y, trans_x_angle_idx, trans_y_angle_idx, phase_trans_x, phase_trans_y, phase_offset] = getSegmentsPhases(x,y,xwt,ywt,Transitions,nTrans,phase_x,phase_y)
    
       angle_x = round(radtodeg(x));
       angle_y = round(radtodeg(y));
       
       all_trans = Transitions(find(x(Transitions)<0));
    
       rand_trans = all_trans(randperm(length(all_trans)));
       
       new_trans = rand_trans(1:nTrans);
 
       for i=1:length(new_trans)
          
           phase_trans_x(i) = angle_x(new_trans(i));
           phase_trans_y(i) = angle_y(new_trans(i));
           
           phase_value_x = phase_trans_x(i) - phase_x;
           
           phase_value_y = phase_trans_y(i) - phase_y;
           
           if phase_x == 0
               
               phase_idx_x = 0;
           
           else
               
               angle_x_tmp = angle_x;
               angle_x_tmp((new_trans(i)+1):end) = [];
               phase_idx_x = max(find(angle_x_tmp==phase_value_x));

               if isempty(phase_idx_x)

                   min_value = min(angle_x_tmp-phase_value_x);
                   phase_idx_x = max(find((angle_x_tmp-phase_value_x)==min_value));

               end
               
           end
           
           if phase_y == 0
               
               phase_idx_y = 0;
               
           else
               
                angle_y_tmp = angle_y;
                angle_y_tmp((new_trans(i)+1):end) = [];
                phase_idx_y = max(find(angle_y_tmp==phase_value_y));
           
                if isempty(phase_idx_y)
               
                    min_value = min(angle_y_tmp-phase_value_y);
                    phase_idx_y = max(find((angle_y_tmp-phase_value_y)==min_value));
               
                end
               
           end
          
%            y_jumps_idx = find(ywt<-6);
%            
%            diff_trans = new_trans(i) - y_jumps_idx;
%            diff_trans(diff_trans>0) = [];
%           
%            min_jump_y_idx = abs(max(diff_trans)) + new_trans(i);
%            
%            x_jumps_idx = find(xwt<-6);
%            
%            diff_trans = new_trans(i) - x_jumps_idx;
%            diff_trans(diff_trans>0) = [];
%           
%            min_jump_x_idx = abs(max(diff_trans)) + new_trans(i);
%            
%            end_segment = min(min_jump_x_idx,min_jump_y_idx);
            
            angle_x_increase = cumsum(diff(x)+2*pi);
            angle_y_increase = cumsum(diff(y)+2*pi);

%             angle_x_increase = angle_x;
%             angle_y_increase = angle_y;
          
           if (phase_idx_y~=0) && (phase_idx_x~=0)
         
               phase_offset = phase_idx_y - phase_idx_x;
               
               end_segment = new_trans(i) + (new_trans(i) - min(phase_idx_x,phase_idx_y));
               
           elseif phase_idx_y == 0
               
               phase_idx_y = phase_idx_x;
               
               phase_offset = 0;
               
               end_segment = new_trans(i) + (new_trans(i) - phase_idx_y);
               
           elseif phase_idx_x == 0
               
               phase_idx_x = phase_idx_y;
               
               phase_offset = 0;
               
               end_segment = new_trans(i) + (new_trans(i) - phase_idx_x);
               
           end
           
           start_x = phase_idx_x;
           start_y = phase_idx_y;
               
           y_vector = start_y:end_segment;
           x_vector = start_x:end_segment-phase_offset;
             
           seg_x(i).phases = angle_x_increase(x_vector);
           seg_y(i).phases = angle_y_increase(y_vector);
           
           phase_increase_trans_x = angle_x_increase(new_trans(i));
           phase_increase_trans_y = angle_y_increase(new_trans(i));
           
%            seg_x(i).phases = seg_x(i).phases - phase_increase_trans_x;
%            seg_y(i).phases = seg_y(i).phases - phase_increase_trans_y;

%             seg_x(i).phases = seg_x(i).phases./10.^floor(log10(seg_x(i).phases)) - seg_x(i).phases;
%             seg_y(i).phases = seg_y(i).phases./10.^floor(log10(seg_y(i).phases)) - seg_y(i).phases;
            
           trans_x_angle_idx(i) = new_trans(i) - phase_idx_x;
           trans_y_angle_idx(i) = new_trans(i) - phase_idx_y;

       end
       
       %save('getSegmentsPhases.mat');
        
 end

function seg = getSegmentsTrans(vector,Transitions,nTrans,segmentsLength)
        
       Transitions(Transitions<segmentsLength) = [];
    
       rand_trans = Transitions(randperm(length(Transitions)));
       
       new_trans = rand_trans(1:nTrans);
       
       seg = zeros(nTrans,segmentsLength*2+1);
       
       offset = -segmentsLength:segmentsLength;
       
       for i=1:length(new_trans)
          
           seg(i,:) = vector(offset + new_trans(i));
           
       end
        
 end

function [uniques,numUnique] = count_unique(x,option)
%COUNT_UNIQUE  Determines unique values, and counts occurrences
%   [uniques,numUnique] = count_unique(x)
%
%   This function determines unique values of an array, and also counts the
%   number of instances of those values.
%
%   This uses the MATLAB builtin function accumarray, and is faster than
%   MATLAB's unique function for intermediate to large sizes of arrays for integer values.  
%   Unlike 'unique' it cannot be used to determine if rows are unique or 
%   operate on cell arrays.
%
%   If float values are passed, it uses MATLAB's logic builtin unique function to
%   determine unique values, and then to count instances.
%
%   Descriptions of Input Variables:
%   x:  Input vector or matrix, N-D.  Must be a type acceptable to
%       accumarray, numeric, logical, char, scalar, or cell array of
%       strings.
%   option: Acceptable values currently only 'float'.  If 'float' is
%           specified, the input x vector will be treated as containing
%           decimal values, regardless of whether it is a float array type.
%
%   Descriptions of Output Variables:
%   uniques:    sorted unique values
%   numUnique:  number of instances of each unique value
%
%   Example(s):
%   >> [uniques] = count_unique(largeArray);
%   >> [uniques,numUnique] = count_unique(largeArray);
%
%   See also: unique, accumarray

% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2009-03-17

testFloat = false;
if nargin == 2 && strcmpi(option,'float')
    testFloat = true;
end

nOut = nargout;
if testFloat
    if nOut < 2
        [uniques] = float_cell_unique(x,nOut);
    else
        [uniques,numUnique] = float_cell_unique(x,nOut);
    end
else
    try %this will fail if the array is float or cell
        if nOut < 2
            [uniques] = int_log_unique(x,nOut);
        else
            [uniques,numUnique] = int_log_unique(x,nOut);
        end
    catch %default to standard approach
        if nOut < 2
            [uniques] = float_cell_unique(x,nOut);
        else
            [uniques,numUnique] = float_cell_unique(x,nOut);
        end
    end
end

end

function [uniques,numUnique] = int_log_unique(x,nOut)
%First, determine the offset for negative values
minVal = min(x(:));

%Check to see if accumarray is appropriate for this function
maxIndex = max(x(:)) - minVal + 1;
if maxIndex / numel(x) > 1000
    error('Accumarray is inefficient for arrays when ind values are >> than the number of elements')
end

%Now, offset to get the index
index = x(:) - minVal + 1;

%Count the occurrences of each index value
numUnique = accumarray(index,1);

%Get the values which occur more than once
uniqueInd = (1:length(numUnique))';
uniques = uniqueInd(numUnique>0) + minVal - 1;

if nOut == 2
    %Trim the numUnique array
    numUnique = numUnique(numUnique>0);
end
end 

function [uniques,numUnique] = float_cell_unique(x,nOut)

if ~iscell(x)
    %First, sort the input vector
    x = sort(x(:));
    numelX = numel(x);
    
    %Check to see if the array type needs to be converted to double
    currClass = class(x);
    isdouble = strcmp(currClass,'double');
    
    if ~isdouble
        x = double(x);
    end
    
    %Check to see if there are any NaNs or Infs, sort returns these either at
    %the beginning or end of an array
    if isnan(x(1)) || isinf(x(1)) || isnan(x(numelX)) || isinf(x(numelX))
        %Check to see if the array contains nans or infs
        xnan = isnan(x);
        xinf = isinf(x);
        testRep = xnan | xinf;
        
        %Remove all of these from the array
        x = x(~testRep);
    end
    
    %Determine break locations of unique values
    uniqueLocs = [true;diff(x) ~= 0];
else
    isdouble = true; %just to avoid conversion on finish
    
    %Sort the rows of the cell array
    x = sort(x(:));
    
    %Determine unique location values
    uniqueLocs = [true;~strcmp(x(1:end-1),x(2:end)) ~= 0] ;
end

%Determine the unique values
uniques = x(uniqueLocs);

if ~isdouble
    x = feval(currClass,x);
end

%Count the number of duplicate values
if nOut == 2
    numUnique = diff([find(uniqueLocs);length(x)+1]);
end
end


end

