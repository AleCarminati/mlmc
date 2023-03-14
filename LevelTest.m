%% 1D test - Level 0
clc
d = 1;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,true,rfs);
l.Y_l
l = l.updateNumSamples(1);
l.Y_l
l = l.updateNumSamples(1000);
l.Y_l
abs(l.Y_l-mean(l.Y_vec))

%% 1D test - Level > 0
clc
d = 1;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,false,rfs);
l.Y_l
l = l.updateNumSamples(1);
l.Y_l
l = l.updateNumSamples(1000);
l.Y_l

%% 2D test - Level 0
clc
d = 2;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,true,rfs);
l.Y_l
l = l.updateNumSamples(1);
l.Y_l
l = l.updateNumSamples(1000);
l.Y_l

%% 2D test - Level > 0
clc
d = 2;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,false,rfs);
l.Y_l
l = l.updateNumSamples(1);
l.Y_l
l = l.updateNumSamples(1000);
l.Y_l