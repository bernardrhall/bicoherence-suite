function frameCache = loadframecache(framePath)
% LOADFRAMECACHE Load frame file cache information
%
% LOADFRAMECACHE reads the frame file cache information stored in the
% specified file.  The resulting cache structure is used to locate
% frame data during subsequent calls to READFRAMEDATA.
%
% usage: frameCache = loadframecache(framePath)
%
%   framePath     name of frame cache file
%
%   frameCache    frame cache structure
%
% The frame cache file should consist of whitespace delimited ASCII
% text and contains one line for each contiguous data segment with a
% common site, frame type, duration, and directory.  Each line should
% consist of the following six columns.
%
%   * site designator (e.g. 'H' or 'L')
%   * frame file type (e.g. 'RDS_R_L3')
%   * GPS start time of segment
%   * GPS stop time of segment
%   * frame file duration in seconds
%   * full path name of directory
%
% The resulting frame cache structure consists of the following six
% fields.
%
%   .sites        cell array of segment site designators
%   .frameTypes   cell array of segment frame file types
%   .startTimes   vector of segment GPS start times
%   .stopTimes    vector of segment GPS stop times
%   .durations    vector of segment frame file durations
%   .directories  cell array of segment directory names
%
% The data segments are inclusive of the specified start time, but
% exclusive of the specified stop time, such that the segment duration
% is simply the difference between the stop and start times.
%
% See also READFRAMEDATA, CREATEFRAMECACHE.pl, and CONVERTLALACHE.pl.

% Shourov K. Chatterji
% shourov@ligo.mit.edu

% $Id: loadframecache.m,v 1.2 2007/08/01 21:06:00 shourov Exp $

% verify correct number of input arguments
error(nargchk(1, 1, nargin));

% read requested frame cache file
[frameCache.sites, frameCache.frameTypes, frameCache.startTimes, ...
 frameCache.stopTimes, frameCache.durations, frameCache.directories] = ...
    dataread('file', framePath, '%s %s %u %u %u %s');
