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

