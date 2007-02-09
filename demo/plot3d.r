require(sna)
data(fauxhigh)
fhd=component.dist(fauxhigh)
sel=which(fhd$mem==1)
fh2=fauxhigh[sel,sel]
m=outer(fauxhigh %v% "Sex" , fauxhigh %v% "Sex", "==")
m=m+1
gplot3d(fh2,vertex.col=(fauxhigh %v%"Grade")[sel],displaylabels=T,
label=(fauxhigh %v%"Grade")[sel],label.col=1,vertex.rad=2,edge.col=m)
   wait=function(a){
     x=proc.time()[3]
     while(proc.time()[3]<=x+a){}
   }
while(T)for(i in 1:360) { rgl.viewpoint(i,abs(180-i)-90,zoom=1.5);
wait(.04) }
