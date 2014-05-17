function [ X,Y ] = getImageNetXY( )

NUM_SAMPLES =100;
NUM_IMAGES = 10;
NUM_FEATURES = 1000;

X=zeros(NUM_SAMPLES*NUM_IMAGES,NUM_FEATURES);
Y=zeros(NUM_SAMPLES*NUM_IMAGES, 1);
counter = 0;

load('wheelchair-n04576002.sbow.mat');
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('television-n04404412.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('persian-cat-n02123394.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('park-bench-n03891251.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('motor-scootern03791053.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('goblet-n03443371.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('french-horn-n03394916.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('fire-engine-n03345487.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('elephant-02504458.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

load('cello-n02992211.sbow.mat');
counter=counter+1
sample = randsample(length(image_sbow),NUM_SAMPLES);
images = image_sbow(sample)
X(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES,:) = getX(images, NUM_SAMPLES, NUM_FEATURES);
Y(counter*NUM_SAMPLES+1:(counter+1)*NUM_SAMPLES)=counter+1;

end

