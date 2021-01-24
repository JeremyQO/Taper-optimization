function ta=TA(p,L0,vf,vb)

t0=(L0/vf)*log((2*vb+vf)/(2*vb-vf));

ta=t0-(L0/vf)*log((vb+0.5*vf-p*vf/L0)/(vb-0.5*vf));