%% Simple prediction of COVID-19 cases
% This script downloads the NY Times COVID-19 data an fits it assuming the
% daily cases are normally distributed. I am no statistician, so this might
% be a really bad fitting algorithm. Take it with a grain of salt.

% Choose location
state='Michigan';

%% Gather data
% Download files    
outfilename = websave('us-states.csv','https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv');
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
tbl = readtable("/home/pnariyoshi/Downloads/us-states.csv", opts);
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
% Get total cases
c=cases(fips==nstate); %c=c(1:end-1);
% Conver to daily cases
dc=diff([0; c]);
% Smooth out a bit
dc=smooth(dc,3);
% Apply log
lc=log(dc);

% How many days to count;
% Ignore the first nd0 days;
nd0=15; % For Michigan 15 gives the day the stay-at-home order started
y=lc(nd0+1:end);
n=-numel(y):-1;

% Create quadratic regression coefficients
x=[ones(numel(y),1) n' n'.^2];

% Estimate parameters
a=x\y;

% Predict next 9 days
nn=(-1:0.01:9)';
xx=[ones(size(nn)) nn nn.^2];
figure(1)
title('Log cases per day')
plot(n,y,n,x*a,xx(:,2),xx*a)

figure(2)
plot(n,exp(y),n,exp(x*a),xx(:,2),exp(xx*a))
title('Cases per day')
legend('Cases','Fit','Prediction','Location','best');
n0=-numel(c):-1;
dnn=0.01;
nnn=(-180:dnn:180)';
xxx=[ones(size(nnn)) nnn nnn.^2];

figure(3)
plot(n0,c,nnn,cumsum(exp(xxx*a))*dnn,nnn,0.9*sum(exp(xxx*a))*dnn*ones(size(nnn)),nnn,0.95*sum(exp(xxx*a))*dnn*ones(size(nnn)),nnn,0.99*sum(exp(xxx*a))*dnn*ones(size(nnn)));
xlim([-numel(c) 14])
legend('Cases','Prediction','90% line','95% line','99% line','Location','best');
grid on; grid minor

% Today
today=sum(exp(xxx(nnn<0,:)*a)*dnn)

% How many next day
tomorrow=sum(exp(xxx(nnn<1,:)*a)*dnn)

% How many next day
inaweek=sum(exp(xxx(nnn<7,:)*a)*dnn)

% Final estimate
final=sum(exp(xxx*a)*dnn)
