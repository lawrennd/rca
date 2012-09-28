function [voteMatrices, Y] = importEVfile(EVfileToRead, allCountryNames, direction)
% function [year, givers, receivers, points] = importEUfile(EVfileToRead)

% IMPORTEUFULE Imports data from the specified Eurovision votes file.
%
% FORMAT
% DESC Converts data from the Eurovision file into a numeric matrix of the
% votes that each country gave to the other countries in a particular year
% of the competition.
%
% ARG EVfileToRead : file to read.
%
% ARG allCountryNames : a cell array of all country name strings,
% aphabetically ordered.
%
% RETURN voteMatrices : a cell array of matrices, whose length is the
% number of years of competition. Each matrix V = voteMatrices{k} is a
% square non-symmetric vote-matrix, where v_ij is the number of points
% awarded from country i to country j, during year k.
%
% RETURN direction : a character string. 'active' returns a numeric matrix
% Y of row size [number of countries]x[years of competition] and column
% size [number of countries]. Each row constists of the votes that a
% particular country gave to every other country (from a complete
% alphatetically ordered list of countries) in a particular year of the
% competition, and it has the following format:
%
%   ([votes to Albania], [votes to Andorra], ..., [votes to United Kingdom])
%
% Every country always rewards the maximum allowed points to itself. This
% encodes for the sensible assumption that a country always likes its own
% song, and the whole row forms an affinity vector of the country towards
% all countries (including itself) for the duration of one competition.
%
% 'passive' returns Y of the same size, but each each row consists of the
% votes that a particular country received from every country (including
% maximum points from itself). If direction is unset then no Y is returned.
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

if nargin > 2 && ~strcmpi(direction, 'active') && ~strcmpi(direction, 'passive')
    error('direction must be ''active'' or ''passive''.')
else
    Y = [];
end

% Import file.
newData = importdata(EVfileToRead, ',', 1);

% Giver and receiver countries are converted to IDs.
years = newData.textdata(2:end,1);  years = cell2mat(years);    years = str2num(years(:,1:4)); %#ok<ST2NM>
givers = countryToID(allCountryNames, newData.textdata(2:end,3));
receivers = countryToID(allCountryNames, newData.textdata(2:end,2));
points = newData.data;
uniqueYears = unique(years);
voteMatrices = cell(length(uniqueYears),1);

for i = 1:length(uniqueYears)
    voteMatrices{i} = 12 * eye(length(allCountryNames));                    % Initialise vote-matrix of this year.
    yearIndx = (years == uniqueYears(i));
    giversThisYear = givers(yearIndx);
    receiversThisYear = receivers(yearIndx);
    pointsThisYear = points(yearIndx);
    uniqueGiversThisYear = unique(giversThisYear);
    for j = 1:length(uniqueGiversThisYear)
        giverThisYearIndx = (giversThisYear == uniqueGiversThisYear(j));
        receiversThisYearThisGiver = receiversThisYear(giverThisYearIndx);
        pointsThisYearThisGiver = pointsThisYear(giverThisYearIndx);
        voteMatrices{i}(uniqueGiversThisYear(j), receiversThisYearThisGiver) = pointsThisYearThisGiver;
    end
end

if nargin > 2
    if strcmpi(direction, 'active')
        for i = 1:length(uniqueYears)
            Y = [Y; voteMatrices{i}]; %#ok<*AGROW>
        end
        Y = Y(sum(Y>0, 2) > 1, :);                                          % Remove empty rows (non-participating countries).
    elseif strcmpi(direction, 'passive')
        for i = 1:length(uniqueYears)
            Y = [Y; voteMatrices{i}'];
        end
%         Y = Y(sum(Y>0) > 1, :); % Remove empty columns (non-participating countries).
    end
end


function countryID = countryToID(allCountryNames, countryName)

% COUNTRYTOID Maps country names to IDs.
%
% FORMAT
% DESC Maps a country name string to a unique ID number. The ID number is
% defined as the index of countryName in the allCountryNames cell
% array. Since the list of countries is fixed and alphabetically ordered,
% the ID of each country is always unique. If countryName is not in
% allCountryNames, zero is returned. If countryName is not just a
% character string, but instead a cell array of country name strings, this
% operation is repeated across all elements in countryName.
%
% ARG allCountryNames : a cell array of all country name strings, aphabetically ordered.
%
% ARG countryName : a country name string (or cell array thereof).
%
% RETURN countryID : a numeric vector of IDs.
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2012
%
% RCA

% Make sure countries are alphabetically ordered.
allCountryNames = unique(allCountryNames);

% Convert singleton to a cell-list.
if ischar(countryName)
    countryName = {countryName};
end

% Mappings between countryName and uniqueCountries and vice versa.
[uniqueCountries, wholeToUniqueInx, uniqueToWholeIndx] = unique(countryName, 'first'); %#ok<ASGLU>

countryID = zeros(length(countryName), 1);
% For each unique country name, find its ID and fill in the corresponding
% indices in countryID.
for uniqIdx = 1:length(uniqueCountries)
    ID = find(cellfun(@(x)(strcmpi(x,uniqueCountries{uniqIdx})), allCountryNames), 1);
    if ~isempty(ID)
        countryID(uniqueToWholeIndx == uniqIdx) = ID;
    end
end



