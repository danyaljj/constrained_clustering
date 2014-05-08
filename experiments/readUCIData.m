function [X, Y] = readUCIData( strClass )
X = [];
Y = [];

if strcmp(strClass, 'balance')
    data = textread( ['../../Data/UCI/balance/balance-scale.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; tmp{1}{1}];
        X = [X; str2num( tmp{1}{2} ) str2num( tmp{1}{3} ) str2num( tmp{1}{4} ) str2num( tmp{1}{5} )  ];
    end
elseif strcmp(strClass, 'glass')
    data = textread( ['../../Data/UCI/glass/glass.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; tmp{1}{11}];
        Xtmp = [];
        for j = 1:10
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'iris')
    data = textread( ['../../Data/UCI/iris/iris.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)-1
        tmp = regexp(data(i), ',', 'split');
        if strcmp( 'Iris-setosa', tmp{1}{5}) 
            Y = [Y; 1];
        elseif strcmp( 'Iris-versicolor', tmp{1}{5})
            Y = [Y; 2];
        elseif strcmp( 'Iris-virginica', tmp{1}{5}) 
            Y = [Y; 3];
        else 
            error('class undefined ! ');
        end 
        Xtmp = [];
        for j = 1:4
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'vertebral')
%     data = textread( ['../../Data/UCI/vertebral/column_2C_weka.arff' ] , '%s', 'delimiter', '\n');
%     data = textread( ['../../Data/UCI/vertebral/column_3C_weka.arff' ] , '%s', 'delimiter', '\n');
    data = textread( ['../../Data/UCI/vertebral/column_3C.dat' ] , '%s', 'delimiter', '\n');
%     data = textread( ['../../Data/UCI/vertebral/column_2C.dat' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ' ', 'split');
        if strcmp( 'DH', tmp{1}{7}) 
            Y = [Y; 1];
        elseif strcmp( 'SL', tmp{1}{7})
            Y = [Y; 2];
        elseif strcmp( 'NO', tmp{1}{7}) 
            Y = [Y; 3];
        else 
            error('class undefined ! ');
        end 
        Xtmp = [];
        for j = 1:6
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
    
elseif strcmp(strClass, 'wdbc')
    
    data = textread( ['../../Data/UCI/wdbc/wdbc.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        if strcmp( 'M', tmp{1}{2}) 
            Y = [Y; 1];
        elseif strcmp( 'B', tmp{1}{2})
            Y = [Y; 2];
        else 
            error('class undefined ! ');
        end 
        Xtmp = [];
        for j = 3:32
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
    
elseif strcmp(strClass, 'wine')
    
    data = textread( ['../../Data/UCI/wine/wine.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; str2num(tmp{1}{1}) ];
        Xtmp = [];
        for j = 2:14
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
end

end