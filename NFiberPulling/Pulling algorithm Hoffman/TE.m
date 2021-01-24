function te=TE(p,L0,vf,vb,Ln)

t0=(L0/vf)*log((2*vb+vf)/(2*vb-vf));
zan=L0+Ln-Ln*0.5*(vf/vb)-t0*(vb-0.5*vf);

%te=t0-(p-Ln)/(vb-0.5*vf);

te=t0-(p-zan)/(vb-0.5*vf);


