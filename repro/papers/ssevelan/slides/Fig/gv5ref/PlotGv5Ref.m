

%% Plot 2D Data;

xrec = 0:0.01:3;
xtk = 0:1:3; ytk=0:0.5:1.5; fontsize=15;
xi = x/1000; zi = z/1000;

%% Plot velocity    
Position = [500, 300,560, 290]; modscale = [1.4 2.5];
Plot2DImageModel(xi,zi,reshape(v*1e-3,n),xtk,ytk,fontsize,Position, modscale);
% hold on;
% plot(model.xsrc/1000,model.zsrc/1000,'g*','MarkerSize',10)
% plot(model.xrec/1000,model.zrec/1000,'bv','MarkerSize',5)


%% Plot Data 
xtkd = [0.5 1 1.5 2 2.5];ytkd = [0 0.5 1 1.5];
TdPos = [500, 200,400, 400]; tdscale = [-30 40]; xlabelname = 'Receiver (km)';
subidx = 1:451;
for pp = 1:nsrc
    rcv = nd(1)*(pp-1)+[1:nd(1)];
    Plot2DTimeGather(xrec,t(subidx),-Dpt(rcv,subidx)',xtkd,ytkd,fontsize,TdPos,tdscale,xlabelname)
    Plot2DTimeGather(xrec,t(subidx),-Sft(rcv,subidx)',xtkd,ytkd,fontsize,TdPos,tdscale,xlabelname)
end

xtkf = [0.5 1 1.5 2 2.5];ytkf = [-1.5 -1 -0.5 0 0.5 1 1.5];

delta_t  = t-2;
DeltaSub = find(abs(delta_t)<=1.5);
FdPos = [500, 200,400, 400]; fdscale = [-30 40]*1e-2; xlabelname = 'Distance (km)';
 Plot2DTimeGather(xrec,delta_t (DeltaSub),F0t(:,DeltaSub)',xtkf,ytkf,fontsize,FdPos,fdscale,xlabelname)