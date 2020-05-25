function cleanResults =...
    cleanEyeTraces(results)
% written by Tiago Marques

%% Inputs
max_mov=5;
pupil_range=[30 100];
sig_mult=3;
mov_avg=60;
n_filt=121;
k_filt=5;
saccade_min=3;
peak_distance=60;
n_rem=1;
max_pupil_disp=25;
radius_cr_med=15;
max_diff_radius_cr=3;

%% Pupil Data
radius_pupil=double(results.pupil_radius);
x_pupil=double(results.pupil_position(:,1));
y_pupil=double(results.pupil_position(:,2));
x_pupil_med=nanmedian(x_pupil);
y_pupil_med=nanmedian(y_pupil);

%% CR Radius
radius_cr=double(results.cr_radius);
%radius_cr_med=median(radius_cr);
radius_cr_std=nanstd(radius_cr);
x_cr=double(results.cr_position(:,1));
y_cr=double(results.cr_position(:,2));
x_cr_med=nanmedian(x_cr);
y_cr_med=nanmedian(y_cr);
x_cr_std=nanstd(x_cr);
y_cr_std=nanstd(y_cr);

%% Maximum allowed parameters
max_diff_x_cr=sig_mult*x_cr_std;
max_diff_y_cr=sig_mult*y_cr_std;

eyeSamples=length(radius_pupil);

%% Initializes error flags
error_flag=false(eyeSamples,1);
error_cr=false(eyeSamples,1);
error_radius=false(eyeSamples,1);
error_range=false(eyeSamples,1);
error_tracker_data=false(eyeSamples,1);
error_pupil_pos=false(eyeSamples,1);
last_good_radius(1)=1;

for i=2:eyeSamples
    diff_pupil_radius=abs(radius_pupil(i)-radius_pupil(last_good_radius(i-1)));
    max_mov=max_mov*(i-last_good_radius(i-1))^(1/3);
    diff_cr_radius=abs(radius_cr(i)-radius_cr_med);
    diff_cr_x= abs(x_cr(i)-x_cr_med); 
    diff_cr_y= abs(y_cr(i)-y_cr_med);
    diff_pupil_y= abs(y_pupil(i)-y_pupil_med);
    diff_pupil_x=abs(x_pupil(i)-x_pupil_med);
    if (diff_cr_radius>max_diff_radius_cr)||(diff_cr_x>max_diff_x_cr)||(diff_cr_y>max_diff_y_cr);
        error_cr(i)=true;
        error_flag(i)=true;
    end
    if (diff_pupil_radius>max_mov)
        error_radius(i)=true;
        error_flag(i)=true;
    end
    if (radius_pupil(i)<pupil_range(1))||(radius_pupil(i)>pupil_range(2))
        error_range(i)=true;
        error_flag(i)=true;
    end
    if isnan(radius_pupil(i))
        error_tracker_data(i)=true;
        error_flag(i)=true;
    end
    if (diff_pupil_y>max_pupil_disp)||(diff_pupil_x>max_pupil_disp)
        error_flag(i)=true;
        error_pupil_pos(i)=true;
    end

    if (i>n_rem)&&(i<(eyeSamples-n_rem))
        if sum(error_flag((i-n_rem):i))==0
            last_good_radius(i,1)=i;
        elseif i>n_rem+1
            last_good_radius((i-n_rem):i,1)=last_good_radius(i-n_rem-1);
        else
            last_good_radius((i-n_rem):i,1)=last_good_radius(1);
        end
    elseif (i<=n_rem)
        if sum(error_flag(1:i))==0
            last_good_radius(i,1)=i;
        else
            last_good_radius(1:i,1)=last_good_radius(1);
        end
    else
        if sum(error_flag((i-n_rem):i))==0
            last_good_radius(i,1)=i;
        else
            last_good_radius((i-n_rem):i,1)=last_good_radius(i-n_rem-1);
        end
    end
end

%% Removes data points next to errors
error_flag_nb=error_flag;
for i=1:n_rem
    error_flag_neg=circshift(error_flag,-i);
    error_flag_neg((end-i+1):end)=0;
    error_flag_pos=circshift(error_flag,i);
    error_flag_pos(1:i)=0;
    error_flag_nb=error_flag_nb|error_flag_neg|error_flag_pos;
end

radius_pupil(error_flag_nb)=NaN;
x_pupil(error_flag_nb)=NaN;
y_pupil(error_flag_nb)=NaN;
x_cr(error_flag_nb)=NaN;
y_cr(error_flag_nb)=NaN;

radius_pupil_fix=fixgaps_alt(radius_pupil);
x_pupil_fix=fixgaps_alt(x_pupil);
y_pupil_fix=fixgaps_alt(y_pupil);
x_cr_fix=fixgaps_alt(x_cr);
y_cr_fix=fixgaps_alt(y_cr);

%% Compensates for global eye movements
x_pupil=x_pupil-x_cr;
y_pupil=y_pupil-y_cr;
x_pupil_fix=x_pupil_fix-x_cr_fix;
y_pupil_fix=y_pupil_fix-y_cr_fix;

%% Centers around the median
x_pupil=x_pupil-nanmedian(x_pupil);
y_pupil=y_pupil-nanmedian(y_pupil);

%% Filters the results
radius_pupil_filt=sgolayfilt(radius_pupil_fix,k_filt,n_filt);
x_pupil_filt=sgolayfilt(x_pupil_fix,k_filt,n_filt);
y_pupil_filt=sgolayfilt(y_pupil_fix,k_filt,n_filt);


x_pupil_filt=x_pupil_filt-median(x_pupil_filt);
y_pupil_filt=y_pupil_filt-median(y_pupil_filt);


%% Detects saccades
pos_x_movavg=slidefun('nanmean',mov_avg,x_pupil_filt,'backward');
pos_y_movavg=slidefun('nanmean',mov_avg,y_pupil_filt,'backward');
x_saccades=abs(x_pupil_filt-pos_x_movavg);
y_saccades=abs(y_pupil_filt-pos_y_movavg);
[~,x_saccades]=findpeaks(x_saccades,'MINPEAKHEIGHT',saccade_min,'MINPEAKDISTANCE',peak_distance);
[~,y_saccades]=findpeaks(y_saccades,'MINPEAKHEIGHT',saccade_min,'MINPEAKDISTANCE',peak_distance);


cleanResults=results;
cleanResults.error_flag=error_flag;
cleanResults.error_cr=error_cr;
cleanResults.error_radius=error_radius;
cleanResults.error_range=error_range;
cleanResults.error_tracker_data=error_tracker_data;

cleanResults.clean_pupil_radius=radius_pupil;
cleanResults.clean_pupil_position=[x_pupil,y_pupil];
cleanResults.x_saccades=x_saccades;
cleanResults.y_saccades=y_saccades;

n_trials=length(results.byTrial.nframes);
step=1;
for i=1:n_trials
    current_frames=results.byTrial.nframes{i};
    cleanResults.byTrial.clean_pupil_radius{i}=radius_pupil(step:(step+current_frames-1));
    cleanResults.byTrial.clean_pupil_position{i}(:,1)=x_pupil(step:(step+current_frames-1));
    cleanResults.byTrial.clean_pupil_position{i}(:,2)=x_pupil(step:(step+current_frames-1));
    step=step+current_frames;
end

