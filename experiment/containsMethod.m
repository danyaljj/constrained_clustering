function ret = containsMethod(e, container)
ret = 0;
if size(container, 1) > 0
    assert(size(container, 1) == 1);
    for i=1:size(container, 2)
        ret = strcmp(container{i}, e);
        if ret == 1
            break
        end
    end
end
