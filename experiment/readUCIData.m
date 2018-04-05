function [X, Y] = readUCIData( strClass )
X = [];
Y = [];

if strcmp(strClass, 'balance')
    data = textread( ['../data/UCI/balance/balance-scale.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; tmp{1}{1}];
        X = [X; str2num( tmp{1}{2} ) str2num( tmp{1}{3} ) str2num( tmp{1}{4} ) str2num( tmp{1}{5} )  ];
    end
    Y(find(Y=='B'))='1';
    Y(find(Y=='L'))='2';
    Y(find(Y=='R'))='3';
    Y=str2num(Y);
elseif strcmp(strClass, 'glass')
    data = textread( ['../data/UCI/glass/glass.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; tmp{1}{11}];
        Xtmp = [];
        for j = 1:10
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
    Y=str2num(Y);
elseif strcmp(strClass, 'iris')
    data = textread( ['../data/UCI/iris/iris.data' ] , '%s', 'delimiter', '\n');
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
    %     data = textread( ['../../data/UCI/vertebral/column_2C_weka.arff' ] , '%s', 'delimiter', '\n');
    %     data = textread( ['../../data/UCI/vertebral/column_3C_weka.arff' ] , '%s', 'delimiter', '\n');
    data = textread( ['../data/UCI/vertebral/column_3C.dat' ] , '%s', 'delimiter', '\n');
    %     data = textread( ['../../data/UCI/vertebral/column_2C.dat' ] , '%s', 'delimiter', '\n');
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

    data = textread( ['../data/UCI/wdbc/wdbc.data' ] , '%s', 'delimiter', '\n');
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

    data = textread( ['../data/UCI/wine/wine.data' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; str2num(tmp{1}{1}) ];
        Xtmp = [];
        for j = 2:14
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'ecoli')
    data = textread( ['../data/UCI/ecoli/ecoli.data.txt' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        %data(i)
        aa = regexprep(data(i), '\s+', ' ');
        tmp = regexp(aa, ' ', 'split');
        %tmp

        if strcmp( 'cp', tmp{1}{9})
            Y = [Y; 1];
        elseif strcmp( 'im', tmp{1}{9})
            Y = [Y; 2];
        elseif strcmp( 'imU', tmp{1}{9})
            Y = [Y; 3];
        elseif strcmp( 'om', tmp{1}{9})
            Y = [Y; 4];
        elseif strcmp( 'omL', tmp{1}{9})
            Y = [Y; 5];
        elseif strcmp( 'pp', tmp{1}{9})
            Y = [Y; 6];
        elseif strcmp( 'imL', tmp{1}{9})
            Y = [Y; 7];
        elseif strcmp( 'imS', tmp{1}{9})
            Y = [Y; 8];
        else
            error('class undefined ! ');
        end

        Xtmp = [];
        for j = 2:8
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'prima')
    data = textread( ['../data/UCI/prima/pima-indians-diabetes.data.txt' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; str2num(tmp{1}{9}) ];
        Xtmp = [];
        for j = 1:8
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'isomosphere')
    data = textread( ['../data/UCI/ionosphere/ionosphere.data.txt' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        if strcmp( 'g', tmp{1}{35})
            Y = [Y; 1];
        elseif strcmp( 'b', tmp{1}{35})
            Y = [Y; 2];
        else
            error('class undefined ! ');
        end

        Xtmp = [];
        for j = 1:34
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'pendigits')
    data = textread( ['../data/UCI/pendigits/pendigits.tra' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; str2num(tmp{1}{17}) ];
        Xtmp = [];
        for j = 1:16
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
    data = textread( ['../data/UCI/pendigits/pendigits.tes' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; str2num(tmp{1}{17}) ];
        Xtmp = [];
        for j = 1:16
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
elseif strcmp(strClass, 'letter.recognition')
    data = textread( ['../data/UCI/letter.recognition/letter-recognition.data.txt' ] , '%s', 'delimiter', '\n');
    for i = 1:size(data,1)
        tmp = regexp(data(i), ',', 'split');
        Y = [Y; tmp{1}{1}-64 ];
        Xtmp = [];
        for j = 2:17
            Xtmp = [Xtmp str2num( tmp{1}{j} ) ];
        end
        X = [X;  Xtmp  ];
    end
end
