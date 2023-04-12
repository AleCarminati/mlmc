%% 1D test - Level 0
clc
d = 1;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,true,rfs);
mean(l.Y_vec)
l = l.updateNumSamples(1);
mean(l.Y_vec)
l = l.updateNumSamples(1000);
mean(l.Y_vec)

%% 1D test - Level > 0
clc
d = 1;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,false,rfs);
mean(l.Y_vec)
l = l.updateNumSamples(1);
mean(l.Y_vec)
l = l.updateNumSamples(1000);
mean(l.Y_vec)

%% 2D test - Level 0
clc
d = 2;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,true,rfs);
mean(l.Y_vec)
l = l.updateNumSamples(1);
mean(l.Y_vec)
l = l.updateNumSamples(1000);
mean(l.Y_vec)

%% 2D test - Level > 0
clc
d = 2;
rfs = RandomFieldSampler(500,1,0.3,d);
l = Level(100,d,16,false,rfs);
mean(l.Y_vec)
l = l.updateNumSamples(1);
mean(l.Y_vec)
l = l.updateNumSamples(1000);
mean(l.Y_vec)