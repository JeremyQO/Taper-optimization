%%%%%%%%%%%%%%%%%%% PULLING ALGORITH DRUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%%%%%%%%%%%%%%%%%%% SECTION 1: CALCULATE PULL PROFILE %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% SUBSECTION A: DEFINITION OF THE TARGET PROFILE %%%%%%%%%%%%%%%%%%

%From r0 to r1: linear taper (angle1)
%From r1 to r_exp: linear taper (angle2)
%From r_exp to rw: exponential taper

%For a single angle taper (angle1=angle2), r1 can be chosen arbitrarily
%between r0 and r_exp

r0=0.0625; %initial radius in mm
L0=0.62; %size of the flame in mm
vb=2.0; %velocit y of the burner in mm/s
Lw=10; %length of the waist in mm
rw=0.001; %waist radius in mm
omega=[0.0001,0.0001]; %taper linear angles in rad [angle2 angle1]
r_exp=0.004; %initial radius of the exponential
r1=0.033; %transition radius between the two linear parts


%%%%%%%%%%%%%% SUBSECTION B: ALGORITHM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=10000; % Number of points for the simulation

%maximum acceleration for the ramps (for XML 210 stages from Newport)
payload = 20;%in kg load on top of motors overestimated to be safe
Mcar = 20;%in kg
amax = 97232/(10*(payload+Mcar));%Max acceleration of the motors
amax = 200;

%%%%%%%%%%%%%%  SUBSECTION C:INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%during the calculation, we need to follow the evolution of the following
%parameters:
cn=[];%radius ratios between two consecutive steps
vf=[];%pulling velocity (mm/s)
rwaist=[];%radius of the waist
Lwn=[];%length of the waist
Ln=[];%traveling length of the flame
t0=[];%heating time of a point totally swept by the flame
zt0=[];%traveling length of a point totally swept by the flame
zan=[];%position of the last point totally swept by the flame
zend=[];%position of the end of the fiber
zta=[];%traveling length of a point that ends in the flame
Lue=[];%difference between extension length of the taper and L0

rwaist(1)=rw; Lwn(1)=Lw;%we initialize to the final dimensions  of the
%waist and reconstruct the fiber backwards 


%We want a continuous transition between the linear part and the
%exponential part of the taper. In order to match the slopes, we fix the
%relative variation of the radius between two consecutive pulling steps for
%the exponential.
c_test=1-sqrt(-L0*omega(1)/r_exp+sqrt(1+(L0*omega(1)/r_exp)^2));
steps_test=(log(rw/r_exp)/log(1-c_test));
steps_exp=floor(steps_test); %the number of pulling steps has to be rounded
%to the closest integer
c_exp=1-(rw/r_exp)^(1/steps_exp);%relative variation of the radius between
%two consecutive steps after rounding


%%%%%%%%%%%%%%%  SUBSECTION D:EXPONENTIAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here, we build the exponential part. As seen below, the waist radii ratios
%between two steps is constant. Consequently, vf, to and zt0 are also
%constant

vf_exp=2*vb*(1-(1-c_exp)^2)/(1+(1-c_exp)^2);
t0_exp=(L0/vf_exp)*log((2*vb+vf_exp)/(2*vb-vf_exp));
zt0_exp=vb*t0_exp-L0;

Lue(1)=-1; %Not used in the exponential part
cn(1)=c_exp;
vf(1)=vf_exp;
t0(1)=t0_exp;
zt0(1)=zt0_exp;
Ln(1)=(Lwn(1)+(vb+0.5*vf(1))*t0(1))/(1+0.5*vf(1)/vb);
zan(1)=L0+Ln(1)-Ln(1)*0.5*(vf(1)/vb)-t0(1)*(vb-0.5*vf(1));
zend(1)=L0+Ln(1)-Ln(1)*0.5*(vf(1)/vb);
zta(1)=-0.5*vf(1)*(Ln(1)/vb)+(vb+0.5*vf(1))*t0(1);

for w=2:steps_exp+1;
    Lue(w)=-1;
    cn(w)=c_exp;
    vf(w)=vf_exp;
    t0(w)=t0_exp;
    zt0(w)=zt0_exp;
    rwaist(w)=rwaist(w-1)./(1-c_exp);
    Lwn(w)=((rwaist(w-1)/rwaist(w))^2)*Lwn(w-1) +(vb-vf_exp/2)*t0_exp;
    Ln(w)=(Lwn(w)+(vb+0.5*vf(w))*t0(w))/(1+0.5*vf(w)/vb);
    zan(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb)-t0(w)*(vb-0.5*vf(w));
    zend(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb);
    zta(w)=-0.5*vf(w)*(Ln(w)/vb)+(vb+0.5*vf(w))*t0(w);
end


%%%%%%%%%%%%%%% SUBSECTION E: LINEAR TAPER ZONE 1 - ANGLE 2 %%%%%%%%%%%%%%%%%%%%

%Now we start building the first part of the linear taper (angle2).
%We first need to reinitialize the parameters at steps_exp+1

rwaist(steps_exp+1)=r_exp;
Lwn(steps_exp+1)=Ln(steps_exp)-vf(steps_exp)*Ln(steps_exp)/(2*vb);
Lue(steps_exp+1)=L0;

rwaist(steps_exp+2)=Lue(steps_exp+1)*tan(omega(1)/2)+rwaist(steps_exp+1);
%Here, we take omega/2 because the taper has already be started in the
%previously.

cn(steps_exp+1)=(rwaist(steps_exp+2)-rwaist(steps_exp+1))/rwaist(steps_exp+2);
%the radii ratios are not constant anymore. We need to redifine it at
%steps_exp+1. Consequently, we need to recalculate the following
%parameters:
vf(steps_exp+1)=2*vb*(1-(1-cn(steps_exp+1))^2)/(1+(1-cn(steps_exp+1))^2);
t0(steps_exp+1)=(L0/vf(steps_exp+1))*log((2*vb+vf(steps_exp+1))/(2*vb-vf(steps_exp+1)));
zt0(steps_exp+1)=vb*t0(steps_exp+1)-L0;
Ln(steps_exp+1)=(Lwn(steps_exp+1)+(vb+0.5*vf(steps_exp+1))*t0(steps_exp+1))/(1+0.5*vf(steps_exp+1)/vb); 
zan(steps_exp+1)=L0+Ln(steps_exp+1)-Ln(steps_exp+1)*0.5*(vf(steps_exp+1)/vb)-t0(steps_exp+1)*(vb-0.5*vf(steps_exp+1));
zend(steps_exp+1)=L0+Ln(steps_exp+1)-Ln(steps_exp+1)*0.5*(vf(steps_exp+1)/vb);
zta(steps_exp+1)=-0.5*vf(steps_exp+1)*(Ln(steps_exp+1)/vb)+(vb+0.5*vf(steps_exp+1))*t0(steps_exp+1);
Lue(steps_exp+2)=2*L0*log(1-cn(steps_exp+1))/((1-cn(steps_exp+1))^2-1);


w=steps_exp+2;

while rwaist(w)<r1;
    rwaist(w+1)=Lue(w)*tan(omega(1))+rwaist(w-1);
    cn(w)=(rwaist(w+1)-rwaist(w))/rwaist(w+1);
    Lue(w+1)=2*L0*log(1-cn(w))/((1-cn(w))^2-1);
    vf(w)=2*vb*(1-(1-cn(w))^2)/(1+(1-cn(w))^2);
    t0(w)=L0/vf(w)*log((2*vb+vf(w))/(2*vb-vf(w)));
    zt0(w)=vb*t0(w)-L0;
    Lwn(w)=((rwaist(w-1)/rwaist(w))^2)*Lwn(w-1)+(vb-vf(w-1)/2)*t0(w-1);
    Ln(w)=(Lwn(w)+(vb+0.5*vf(w))*t0(w))/(1+0.5*vf(w)/vb);
    zan(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb)-t0(w)*(vb-0.5*vf(w));
    zend(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb);
    zta(w)=-0.5*vf(w)*(Ln(w)/vb)+(vb+0.5*vf(w))*t0(w);
    w=w+1;
end


%%%%%%%%%%%%%%% SUBSECTION F: LINEAR TAPER ZONE 2 - ANGLE 1 %%%%%%%%%%%%%%%%%%%%

%reinitialization of the parameters before starting the new linear section

rwaist(w)=L0*tan(omega(2)/2)+rwaist(w-1);
cn(w-1)=(rwaist(w)-rwaist(w-1))/rwaist(w);
Lue(w)=2*L0*log(1-cn(w-1))/((1-cn(w-1))^2-1);
vf(w-1)=2*vb*(1-(1-cn(w-1))^2)/(1+(1-cn(w-1))^2);
t0(w-1)=L0/vf(w-1)*log((2*vb+vf(w-1))/(2*vb-vf(w-1)));
zt0(w-1)=vb*t0(w-1)-L0;
Lwn(w-1)=((rwaist(w-2)/rwaist(w-1))^2)*Lwn(w-2)+(vb-vf(w-2)/2)*t0(w-2);
Ln(w-1)=(Lwn(w-1)+(vb+0.5*vf(w-1))*t0(w-1))/(1+0.5*vf(w-1)/vb);
zan(w-1)=L0+Ln(w-1)-Ln(w-1)*0.5*(vf(w-1)/vb)-t0(w-1)*(vb-0.5*vf(w-1));
zend(w-1)=L0+Ln(w-1)-Ln(w-1)*0.5*(vf(w-1)/vb);
zta(w-1)=-0.5*vf(w-1)*(Ln(w-1)/vb)+(vb+0.5*vf(w-1))*t0(w-1);

steps1=w-1;

while rwaist(w)<r0;
    rwaist(w+1)=Lue(w)*tan(omega(2))+rwaist(w-1);
    cn(w)=(rwaist(w+1)-rwaist(w))/rwaist(w+1);
    Lue(w+1)=2*L0*log(1-cn(w))/((1-cn(w))^2-1);
    vf(w)=2*vb*(1-(1-cn(w))^2)/(1+(1-cn(w))^2);
    t0(w)=L0/vf(w)*log((2*vb+vf(w))/(2*vb-vf(w)));
    zt0(w)=vb*t0(w)-L0;
    Lwn(w)=((rwaist(w-1)/rwaist(w))^2)*Lwn(w-1)+(vb-vf(w-1)/2)*t0(w-1);
    Ln(w)=(Lwn(w)+(vb+0.5*vf(w))*t0(w))/(1+0.5*vf(w)/vb);
    zan(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb)-t0(w)*(vb-0.5*vf(w));
    zend(w)=L0+Ln(w)-Ln(w)*0.5*(vf(w)/vb);
    zta(w)=-0.5*vf(w)*(Ln(w)/vb)+(vb+0.5*vf(w))*t0(w);
    w=w+1;
end


%%%%%%%%%%%%%%% SUBSECTION G: ENDING THE LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The algorithm ends when the radius becomes bigger than r0.
%We need to correct the last step to stop at r0 exactly.

steps=w-1;
rwaist(steps+1)=r0;
cn(steps)=1-rwaist(steps)/r0;
vf(steps)=2*vb*(1-(1-cn(steps))^2)/(1+(1-cn(steps))^2);
Lue(steps)=2*L0*log(1-cn(steps))/((1-cn(steps))^2-1);
t0(steps)=L0/vf(steps)*log((2*vb+vf(steps))/(2*vb-vf(steps)));
zt0(steps)=vb*t0(steps)-L0;
Lwn(steps)=((rwaist(steps-1)/rwaist(steps))^2)*Lwn(steps-1)+(vb-vf(steps-1)/2)*t0(steps-1);
Ln(steps)=(Lwn(steps)+(vb+0.5*vf(steps))*t0(steps))/(1+0.5*vf(steps)/vb);
zan(steps)=L0+Ln(steps)-Ln(steps)*0.5*(vf(steps)/vb)-t0(steps)*(vb-0.5*vf(steps));
zend(steps)=L0+Ln(steps)-Ln(steps)*0.5*(vf(steps)/vb);
zta(steps)=-0.5*vf(steps)*(Ln(steps)/vb)+(vb+0.5*vf(steps))*t0(steps);


%%%%%%%%%%%%%%%%%%% SECTION 2: SIMULATION OF THE PULL %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% TRAJECTORIES OF THE MOTORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Using the parameters calculated previously, we perform here a full
%simulation of the pull, starting from r0 and ending at rw.

%We calculate the trajectories of the motors and generate the PVTM matrix
%used by the controller that commands the motors
%PVTM=[T x1 v1 x2 v2] where:
%T is the duration of a step
%x1 is the relative position motor1 has to reach at T
%x2 is the relative position motor2 has to reach at T
%v1 is the velocity motor1 has to reach at T
%v2 is the velocity motor2 has to reach at T

cn=fliplr(cn);
vf=fliplr(vf);
t0=fliplr(t0);
zt0=fliplr(zt0);
Lwn=fliplr(Lwn);
rwaist=fliplr(rwaist);
Lue=fliplr(Lue);
Ln=fliplr(Ln);
zan=fliplr(zan);
zend=fliplr(zend);
zta=fliplr(zta);

time = cumsum(Ln./vb);  %defines the time at the end of each step

PVTM=[];


%We generate vectors with each element representing a point of the fiber
%(radius and position)
znew=0:zend(1)/(N-1):zend(1);
rnew=r0+zeros(size(znew));
zflip=znew; %To flip the axis and redefine the zero after each step
zmatrix = zeros(N,steps);
rmatrix = zeros(N,steps);


for j=1:steps
    %Simulation of the pull
    if(j>1)
        zold=-znew+(Ln(j-1)+L0);
    else
        zold=znew;
    end
    rold=rnew;
    
    gate = (zold<0);%points located before the flame
    znew(gate)=zold(gate)-0.5*vf(j)*(Ln(j)/vb);
    rnew(gate)=rold(gate);
    
    gate = (zold>=0) & (zold<L0);%points starting in the flame 
    znew(gate)=zold(gate)-0.5*vf(j)*TNbf(zold(gate),L0,vf(j),vb,Ln(j))+ZTA(TA(zold(gate),L0,vf(j),vb),zold(gate),vb);
    rnew(gate)=rold(gate).*exp(-0.5.*(vf(j)/L0).*TA(zold(gate),L0,vf(j),vb));
    
    gate = (zold>=L0) & (zold<zan(j));%points fully swept by the flame
    znew(gate)=zold(gate)-0.5*vf(j)*TNbf(zold(gate),L0,vf(j),vb,Ln(j))+zt0(j)+0.5*vf(j)*TV(zold(gate),L0,vf(j),vb);
    rnew(gate)=rold(gate).*exp(-0.5*(vf(j)/L0)*t0(j));
    
    gate = (zold>=zan(j)) & (zold<=zend(j));%points ending in the flame
    znew(gate)=L0+Ln(j)+L0*(vb-0.5*vf(j))/vf(j)*(1-exp(vf(j)/L0*(Ln(j)/vb+(zold(gate)-L0)/(0.5*vf(j)-vb))));
    rnew(gate)=rold(gate).*exp(-0.5.*(vf(j)/L0).*(t0(j)+L0/vf(j).*log(1-vf(j)/L0.*(znew(gate)-Ln(j))./(0.5.*vf(j)+vb))));
    
    gate = (zold>zend(j));%points not reached by the flame
    znew(gate)=zold(gate)+0.5*vf(j)*(Ln(j)/vb);
    rnew(gate)=rold(gate);
    
    zInterp = linspace(min(znew(:)),max(znew(:)),N);%Equally spaced z-axis
    rnew = interp1(real(znew),rnew,real(zInterp));
%     rnew = interp1((znew),rnew,(zInterp));
    znew = zInterp;
    
    zflip=znew-min(znew(:)); %Flipping the axis so that the zero in each step is        %redefined to the outer edge of the flame.
    
    zmatrix(:,j)=zflip;
    rmatrix(:,j)=rnew;
    
    
    plot(zflip,rnew,'b');
    ylim([0 0.065]);
    xlabel('Fiber Length [mm]');
    ylabel('Radius [mm]');
    title('Fiber Stretching')
    drawnow;
    %Generation of the PVTM matrix
    t=Ln(j)/vb;%duration of the step
    
    p1=abs(-vf(j)*t/2+vb*t*(-1)^j);%distance motor1 has to travel
    p2=abs(vf(j)*t/2+vb*t*(-1)^j);%distance motor2 has to travel
    vmot1=-vf(j)/2+vb*(-1)^j;%velocity motor1 has to reach
    vmot2=vf(j)/2+vb*(-1)^j;%velocity motor2 has to reach
   
    
    if (mod(j,2)==0) %Even steps. The motors go in one direction
        a=amax;
        dt1 = abs(vmot1/a);
        dx1 = abs(vmot1*vmot1/(2*a));
        dt2 = abs(vmot2/a);
        dx2 = abs(vmot1*(abs(vmot2)-0.5*abs(vmot1))/a);
        vdiff=abs(vmot2)-abs(vmot1);
        
        %At the beiginning and at the end of a step, the velocities of the
        %motors has to be 0. The motors are be ramped up and ramped to a
        %zero-velocity
        %The ramps are designed to be fast and jerkless
        
        ramp_up_1 = [vdiff/a 0 0 0.5*vdiff^2/a vdiff];
        ramp_up_2 = [dt1 -dx1 -vmot1 dx2 vmot2];
        PVT_step = [p2/vmot2 -p1 -vmot1 p2 vmot2];
        ramp_down_1 = [dt1 -dx1 0 dx2 vdiff];
        ramp_down_2 = [vdiff/a 0 0 0.5*vdiff^2/a 0];
        
    else %Odd steps. The motors go in the opposite direction
        a=amax;
        dt1 = abs(vmot1/a);
        dx1 = abs(vmot2*(abs(vmot1)-0.5*abs(vmot2))/a);
        dt2 = abs(vmot2/a);
        dx2 = abs(vmot2*vmot2/(2*a));
        vdiff=abs(vmot1)-abs(vmot2);
        
        ramp_up_1 = [vdiff/a 0.5*vdiff^2/a vdiff 0 0];
        ramp_up_2 = [dt2 dx1 -vmot1 -dx2 vmot2];
        PVT_step = [-p1/vmot1 p1 -vmot1 -p2 vmot2];
        ramp_down_1 = [dt2 dx1 vdiff -dx2 0];
        ramp_down_2 = [vdiff/a 0.5*vdiff^2/a 0 0 0];
        
    end
    
    PVTM(5*j-4,:)=ramp_up_1;
    PVTM(5*j-3,:)=ramp_up_2;
    PVTM(5*j-2,:)=PVT_step;
    PVTM(5*j-1,:)=ramp_down_1;
    PVTM(5*j,:)=ramp_down_2;
    

end

%% Output
PullingTime = sum(PVTM(:,1))

%%


%%%%%%%%%%%%%%% SECTION 3: SAVE FILES OF INTEREST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


saveStr1 = ['PVT_CAT_vb',num2str(vb),'_r0_',num2str(r0*10^3),'_R',num2str(r1*10^3),'_Rexp',num2str(r_exp*10^3),'_rw',num2str(rw*10^3),'um_Lw',num2str(Lw),'mm_Omega',num2str(omega(1)*10^3),'-',num2str(omega(2)*10^3),'mrad.txt'];
saveStr2 = ['rwaist',saveStr1];
saveStr3 = ['time',saveStr1];
saveStr4 = ['rnew',saveStr1];
saveStr5 = ['znew',saveStr1];

dlmwrite(saveStr1,PVTM,'precision','%.6f','newline','pc');
dlmwrite(saveStr2,rwaist.','precision','%.6f','newline','pc');
dlmwrite(saveStr3,time.','precision','%.6f','newline','pc');
dlmwrite(saveStr4,rnew.','precision','%.6f','newline','pc');
dlmwrite(saveStr5,znew.','precision','%.6f','newline','pc');
