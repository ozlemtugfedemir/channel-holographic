%close all
Mh = 64;  %number of antennas along the horizontal axis
Mv = 64;  %number of antennas along the vertical axis
dH = 0.25; %horizontal inter-antenna distance in number of wavelengths
dV = 0.25; %vertical inter-antenna distance in number of wavelengths

rng(111)
N = 20; %number of clusters
rTau = 2.3;
hBS = 25;
hUE = 1.5;
muu = log10(16/3);
d2d = (0.9-muu)/2.1*1000;
elev0 = 90+atand((hBS-hUE)/d2d);
ZODoffset = -10^(-0.62*log10(d2d)+1.93);
DS = 10^(-6.44);
taunPrime = -rTau*DS*log(rand(1,N));
taun = sort(taunPrime-min(taunPrime),'ascend');
PnPrime = exp(-taun*(rTau-1)/(rTau*DS)).*10.^(-3*randn(1,N)/10);
Pn = PnPrime/sum(PnPrime);
ASD = 10^(1.41);
Cazim = 1.289;
varphiPrime = 2*ASD/1.4*sqrt(-log(Pn/max(Pn)))/Cazim;
varphi = sign(randn(1,N)).*varphiPrime + ASD/7*randn(1,N);
ZSD = 2/3*8;
Celev = 1.1918;
thetaPrime = -ZSD*(log(Pn/max(Pn)))/Celev;
theta = sign(randn(1,N)).*thetaPrime + ZSD/7*randn(1,N) + elev0 + ZODoffset;
theta =  90 - theta;
varphi = varphi*pi/180;
theta = theta*pi/180;

%%%%%%%%%


ASD_varphi = deg2rad(2); %angular standard deviation for azimuth angle
ASD_theta = deg2rad(2);  %angular standard deviation for elevation angle
warning('off')

M = Mh*Mv; %number of total antennas
posY = vec(repmat(0:Mh-1,Mv,1)*dH);  %horizontal position of antennas
posZ = vec(repmat((0:Mv-1).',1,Mh)*dV);  %vertical position of antennas

Rexact = zeros(M,M); %initialize exact correlation matrix for non-isotropic case
Rapproximate = zeros(M,M); %initialize approximate correlation matrix for non-isotropic case

%Parameters for cosine antenna pattern from the paper
aExp = 0;
bExp0 = 0;
bExp = bExp0+1;


scalee = 0;
scaleeB = 0;
for n = 1:N
    F = @(delta,epsilon)cos(varphi(n)+delta).^aExp.*cos(theta(n)+epsilon).^bExp.* ...
        exp(cos(2*delta)/(4*ASD_varphi^2)-200).*...
        exp(cos(2*epsilon)/(4*ASD_theta^2)-200);
    scalee = scalee + Pn(n)*integral2(F,-varphi(n)-pi/2,-varphi(n)+pi/2,-theta(n)-pi/2,-theta(n)+pi/2,'RelTol',1e-9);
    
    scaleeB = scaleeB + Pn(n)*(cos(varphi(n))^aExp*cos(theta(n))^bExp);
end


%% Go through all the columns of the first row
for row = 1:M
    row
    
    %Construct the exact spatial correlation matrix using closed-form
    %approximation from the paper
    for column = row+1:M
        
        distanceY = posY(row)-posY(column);
        distanceZ = posZ(row)-posZ(column);
        for n = 1:N
            F = @(delta,epsilon)exp(1i*2*pi*distanceY*sin(varphi(n)+delta).*cos(theta(n)+epsilon)).*...
                exp(1i*2*pi*distanceZ*sin(theta(n)+epsilon)).*...
                cos(varphi(n)+delta).^aExp.*cos(theta(n)+epsilon).^bExp.* ...
                exp(cos(2*delta)/(4*ASD_varphi^2)-200).*...
                exp(cos(2*epsilon)/(4*ASD_theta^2)-200);
            
            
            Rexact(row,column) = Rexact(row,column) + Pn(n)/scalee*...
                integral2(F,-varphi(n)-pi/2,-varphi(n)+pi/2,-theta(n)-pi/2,-theta(n)+pi/2,'RelTol',1e-9);
            
            
            A = exp(1i*2*pi*(distanceY*sin(varphi(n))*cos(theta(n))+distanceZ*sin(theta(n))));
            B = 2*pi*distanceY*cos(varphi(n))*cos(theta(n));
            C = -2*pi*distanceY*cos(varphi(n))*sin(theta(n));
            D = -2*pi*distanceY*sin(varphi(n))*sin(theta(n))+2*pi*distanceZ*cos(theta(n));
            
            sigma2tilde = ASD_varphi^2/(1+C^2*ASD_varphi^2*ASD_theta^2);
            X = cos(varphi(n))^aExp*(cos(theta(n))^bExp-1i*bExp*cos(theta(n))^(bExp-1)*sin(theta(n))*ASD_theta^2*D);
            Y = aExp*cos(varphi(n))^(aExp-1)*sin(varphi(n))*cos(theta(n))^bExp ...
                + 1i*ASD_theta^2*bExp*cos(theta(n))^(bExp-1)*sin(theta(n))*...
                (cos(varphi(n))^aExp*C-aExp*cos(varphi(n))^(aExp-1)*sin(varphi(n))*D);
            Rapproximate(row,column) = Rapproximate(row,column) + Pn(n)*A*sqrt(sigma2tilde)/ASD_varphi*...
                exp(D^2*ASD_theta^2*(C^2*ASD_theta^2*sigma2tilde -1)/2)*...
                (X-Y*(1i*B*sigma2tilde-C*D*ASD_theta^2*sigma2tilde))*...
                exp(-1i*B*C*D*ASD_theta^2*sigma2tilde)*...
                exp(-B^2*sigma2tilde/2)/scaleeB;
        end
        
        
    end
end

%Complete the missing parts of the correlation matrices
Rexact = (Rexact+Rexact')+eye(M);
Rapproximate = (Rapproximate+Rapproximate')+eye(M);


%% Eigenvalue decomposition
[Uexact,Dexact] = eig(Rexact);
[Uapproximate,Dapproximate] = eig(Rapproximate);

%To prevent any numerical issues, equate very small negative eigenvalues
%to zero
Dexact(Dexact<0)= 0;
Dapproximate(Dapproximate<0)= 0;

%Determine the exact rank
dd = sort(diag(Dexact),'descend');
for rr = 1:M
    if sum(dd(1:rr))/sum(dd)>0.99999
        rank_exact = rr;
        break
    end
end


dd = sort(diag(Dapproximate),'descend');
for rr = 1:M
    if sum(dd(1:rr))/sum(dd)>0.99999
        rank_approximate = rr;
        break
    end
end


Rexact = Uexact*Dexact*Uexact';
Rapproximate = Uapproximate*Dapproximate*Uapproximate';

figure
plot(sort(diag(Dexact),'descend'))
hold on
plot(sort(diag(Dapproximate),'descend'))

simil = trace(Rexact'*Rapproximate)/sqrt(trace(Rexact'*Rexact)*trace(Rapproximate'*Rapproximate))

