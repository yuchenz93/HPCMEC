function a = berow(a)

if isempty(a)
    a = a;
    return
end

if ~isvector(a)
   error('input is not a vector') 
end
if isrow(a)
    a = a;
end
if iscolumn(a)
    a = a';
end


end

