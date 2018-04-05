function out = containsMethod(name, titles)
name 
titles 
 
DataSize = length(titles); 
for i=1:DataSize
    if strcmp( titles{i}, name )
        out = true; 
        return; 
    end
end 
out = false; 
end
