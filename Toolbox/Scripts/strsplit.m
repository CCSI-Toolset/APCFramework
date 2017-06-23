function parts = strsplit(str, splitstr)

    str_idx = strfind(str, splitstr);
    splitlen = length(splitstr);
    parts = {};
    i_prev = 1;
    i = 2; flag = 0;
    while i <= length(str)
        if ~isempty(find(i == str_idx, 1)) && ~flag
            parts{end+1} = str(i_prev:i-1);
            i = i + splitlen;
            i_prev = i;
            flag = 1;
        else
            i = i + 1;
            flag = 0;
        end
    end
    if i_prev < i && ~flag
        parts{end+1} = str(i_prev:i-1);
        i_prev = i;
        flag = 1;
    end
