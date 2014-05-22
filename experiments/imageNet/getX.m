function [ X ] = getX( images, NUM_SAMPLES, NUM_FEATURES)

X = [];
eps = .1;

for i=1:1:NUM_SAMPLES
    ints = images(i).sbow.word;
    counts = zeros(1000,1);
    l = length(ints);
    
    for j=1:1:NUM_FEATURES
        ind = find(ints == j);
        counts(j,1) = (length(ind)+eps) / (l+NUM_FEATURES*eps);
    end
    X(end+1,:)=counts;
end


end

