% this file checks whether a given image is poorly copied

function leaveout = ADNIerrors()
    leaveout = cell(3,1);
    leaveout{1} = [06 09 21 31 35 36 47 59];
    leaveout{2} = [002 004 005 006 021 023 052 053 061 068 071 075 ...
                   076 078 093 094 095 098 099 105 126 127 132 133 ...
                   134 136 140 147 153 156 158 164 181 191 197];
    leaveout{3} = [05 55 56 67 81 86];
end

