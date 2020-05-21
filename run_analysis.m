%% Simple prediction of COVID-19 cases
% This script downloads the NY Times COVID-19 data an fits it assuming the
% daily cases are normally distributed. I am no statistician, so this might
% be a really bad fitting algorithm. Take it with a grain of salt.

% Choose location
state='Michigan';

%% Gather data
% Download files    
websave('us-states.csv','https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv');
% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["date1", "state", "fips", "cases", "deaths"];
opts.VariableTypes = ["datetime", "categorical", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, "state", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "date1", "InputFormat", "yyyy-MM-dd");
% Import the data
tbl = readtable("us-states.csv", opts);
% Convert to output type
date1 = tbl.date1;
states = tbl.state;
fips = tbl.fips;
cases = tbl.cases;
deaths = tbl.deaths;
% Clear temporary variables
clear opts tbl

%% Process data
% Get fips number
nstate=fips(find(states==state,1));
% Get total cases and dates
d=date1(fips==nstate);
c=cases(fips==nstate);
c=smooth(c,7);
% Convert to daily cases
dc=diff([0; c]);
% Smooth out a bit
% dc=smooth(dc,7);

% How many days to count;
% Ignore the first nd0 days;
% nd0=1; % First data point
nd0=find(d=='2020-04-08',1); % 
% nd0=find(d=='2020-03-24',1); % For Michigan, gives the day the stay-at-home order started
% ndlast=find(d=='2020-04-12',1); % For Michigan, gives 2 days after the protests started
% nd0=find(d=='2020-04-26',1); % For Michigan, start about 2 weeks after protests started
ndlast=numel(d); % Uses the last available day

ncurrent=numel(d);
y=dc(nd0:ndlast);
n=(nd0-ncurrent:ndlast-ncurrent)';

% Create model
ft=fittype(@(A,mu,sigma,x) A*normpdf(x,mu,sigma), 'independent', 'x', 'dependent', 'y','coefficients',{'A','mu','sigma'});
% Estimate parameters
F0=fit(n,y,ft,'Startpoint',[max(y),-10,20],'Robust', 'on','Upper',[100*max(y),40,100],'Lower',[1,-100,0.1]);
A=F0.A; mu=F0.mu; sigma=F0.sigma;

% Predict next 30 days
nn=(ndlast-ncurrent:0.01:30)';
% Predict -180 days ago to 180 days later
n0=-numel(c)+1:0;
dnn=0.01;
nnn=(-180:dnn:180)';

% Create an offset for the CDF (there are probably better ways to do this
offset = interp1(n0,c,mean(n))-A*normcdf(mean(n),mu,sigma);


figure(1)
subplot(2,4,1)
plot(n,log(y),'x-',n,log(F0(n)),'x-')
title('Log cases per day')
subplot(2,4,5)
plot(n0,log(dc),'x-',n,log(F0(n)),'x-')
title('Log cases per day')

subplot(2,4,2)
plot(n,y,'x-',n,F0(n),'x-',nn,F0(nn))
title('Cases per day')
legend('Cases','Fit','Prediction','Location','best');
subplot(2,4,6)
plot(n0,dc,'x-',n,F0(n),'x-',nn,F0(nn))
title('Log cases per day')
legend('Cases','Fit','Prediction','Location','best');

subplot(1,2,2)
plot(n0,c,nnn,offset+A*normcdf(nnn,mu,sigma),nnn,offset+max(A*normcdf(nnn,mu,sigma))*ones(size(nnn))*0.9,nnn,offset+max(A*normcdf(nnn,mu,sigma))*ones(size(nnn))*0.95,nnn,offset+max(A*normcdf(nnn,mu,sigma))*ones(size(nnn))*0.99);
xlim([-numel(c) 14])
title('Cumulative cases by day')
legend('Cases','Prediction','90% line','95% line','99% line','Location','best');
grid on; grid minor

% Today
today=offset+A*normcdf(0,mu,sigma)

% How many next day
tomorrow=offset+A*normcdf(1,mu,sigma)

% How many next week
inaweek=offset+A*normcdf(7,mu,sigma)

% Final estimate
final=offset+A
