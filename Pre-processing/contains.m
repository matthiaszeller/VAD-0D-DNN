%contains returns true if the pattern is found in the string(s).
%
%  Usage:
%   b = contains(str,pattern)
%   b = contains(str,pattern,'IgnoreCase',ignorecase)
%
%  str and pattern may be a character array or a cell array of strings and do not
%   have to be the same size.
%  b is the same size as str and returns true for each element in str that
%   matches any of the pattern strings specified.
%  if ignorecase is true, then the match is case insensitive (default is false)
%
%   Examples
%       str = 'data.tar.gz';
%       P = 'tar';
%       contains(str,P)                   returns  1
%
%       str = {'abstracts.docx','data.tar.gz'};
%       P = 'tar';
%       contains(str,P)                   returns  [0 1]
%
%       str = 'data.tar.gz';
%       P = {'docx','tar'};
%       contains(str,P)                   returns  1
%
%       str ={'DATA.TAR.GZ','SUMMARY.PPT'};
%       P = 'tar';
%       contains(str,P,'IgnoreCase',true) returns  [1 0]
% SOURCE: https://savannah.gnu.org/bugs/?56065
function [b] = contains(str,pattern,~,ignorecase=false)
    if ignorecase
        str=lower(str);
        pattern=lower(pattern);
    endif
    if ischar(str)
        if ischar(pattern)
            b=~isempty(strfind(str,pattern));
        elseif iscell(pattern)
            b=any(cellfun(@(p)~isempty(p),strfind(str,pattern)));
        else
            error('pattern must be a char array or cell array of strings');
        endif
    elseif iscell(str)
        if isstr(pattern)
            b=cellfun(@(s)~isempty(s),strfind(str,pattern));
        elseif iscell(pattern)
            b=cellfun(@(s)any(cellfun(@(p)~isempty(p),strfind(s,pattern))),str);
        else
            error('pattern must be a char array or cell array of strings');
        endif
    else
        error('str must be a char array or cell array of strings');
    endif
endfunction

%!assert(contains('data.tar.gz','tar'),true)
%!assert(contains('peppers, onions, and mushrooms','onion'),true);
%!assert(contains('peppers, onions, and mushrooms','pineapples'),false);

%!assert(contains({'abstracts.docx','data.tar.gz'},'tar'),logical([0 1]))
%!assert(contains({'abstracts.docx','data.tar.gz'},'zip'),logical([0 0]))

%!assert(contains('data.tar.gz',{'docx','tar'}),true)
%!assert(contains('data.tar.gz',{'7z','zip'}),false)

%!test
%! str = {'Mary Ann Jones','Christopher Matthew Burns','John Paul Smith'};
%! pattern = {'Ann','Paul'};
%! assert(contains(str,pattern),logical([1 0 1]))
%! assert(str(contains(str,pattern)),{'Mary Ann Jones','John Paul Smith'})

%!test
%! str ={'DATA.TAR.GZ','SUMMARY.PPT'};
%! P = 'tar';
%! assert(contains(str,P,'IgnoreCase',true),logical([1 0]))
%! assert(contains(str,P,'IgnoreCase',false),logical([0 0]))