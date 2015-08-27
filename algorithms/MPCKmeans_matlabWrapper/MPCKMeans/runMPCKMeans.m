function [assignment, centroids_vectors] = runMPCKMeans(X, Y, k, E) 

javaaddpath('C:\Users\Daniel\Desktop\Constrained_codes\otherCodes\wekaUT\wekaUT\weka-latest\src\out\artifacts\weka_latest_jar\weka-latest.jar');

import weka.clusterers.MPCKMeans
import weka.clusterers.InstancePair
import weka.core.Instances
import weka.core.Instance
import weka.core.FastVector
import weka.core.Attribute
import java.util.ArrayList

size1 = size(X,1);
featureVector = FastVector();
featureVector.addElement(Attribute('att1'));
featureVector.addElement(Attribute('att2'));
featureVector.addElement(Attribute('label'));

data = Instances('Persons', featureVector, size1);

instance = Instance(3);
for j = 1:size1 % data points 
    for i =0:1 % attributes 
        i 
        attr = featureVector.elementAt(i); 
        instance.setValue(attr, X(j, i+1));
    end
    attr = featureVector.elementAt(2)
    instance.setValue(attr, Y(j));

%     instance.setClassValue(-1^(j));
    data.add(instance)
end
data.setClassIndex(2);

mpckmeans = MPCKMeans();
mpckmeans.setNumClusters(k);

%classIndex = -1;
classIndex = data.numAttributes()-1; 
if classIndex >= 0
    data.setClassIndex(classIndex); %  starts with 0
    clusterData = Instances(data);
    mpckmeans.setNumClusters(k);
    clusterData.deleteClassAttribute();
else
    clusterData = Instances(data);
end
	
% constraints
if 0
    pairs = ArrayList();
    pair = InstancePair(0, 1, InstancePair.CANNOT_LINK);
    pairs.add(pair);
    pair = InstancePair(1, 2, InstancePair.MUST_LINK);
    pairs.add(pair);
    labeledPairs = pairs;
else
    labeledPairs = ArrayList(0);
end 
mpckmeans.setTotalTrainWithLabels(data);
%System.out.println();
numClusters1 = mpckmeans.getNumClusters(); 
assignments1 = mpckmeans.getClusterAssignments(); 

mpckmeans.buildClusterer(labeledPairs, clusterData, data, mpckmeans.getNumClusters(), data.numInstances());
% assignments = zeros(size1, 1); 
nCorrect = 0;
totalTrainWithLabels = mpckmeans.getTotalTrainWithLabels(); 
assignment = mpckmeans.getClusterAssignments();
% for i=1:mpckmeans.getTotalTrainWithLabels.numInstances()
%     assignments(i) = cluster_i; 
% end 

centroids = mpckmeans.getClusterCentroids(); 
centroids_vectors = [];
for i = 0:centroids.numAttributes()-1
    centroids_vectors = [centroids_vectors  centroids.attributeToDoubleArray(i)]; 
end 

end