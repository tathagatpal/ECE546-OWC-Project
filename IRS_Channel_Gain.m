clc; clear all; close all;

rho=0.8;
m=2;
P=20;

%S coordinates
xs=-0.26;ys=2.5;zs=0.5;

%D coordinates
% xd = 2; 
xd = -1:0.1:1;
for ind = 1:length(xd)
yd = 2; zd = 3;
wm=0.1;
hm=0.1;
irradiance_PS_Mi=0;

nm=5;% no of mirrors in each rows and columns
R_kl = {};
S = {}; D = {};

for l=1:nm
    for k=1:nm
        %Center of kth row and lth column  
        R_kl{k,l} = [((xs+(wm/2)+(l-1)*wm));...
            0;...
            ((zs+(hm/2)+(k-1)*hm))]; 
        
        S{k,l} = [-((xs+(wm/2)+(l-1)*wm));...
            ys;...
            -((zs+(hm/2)+(k-1)*hm))];
        
        D{k,l} = [xd(ind)-((xs+(wm/2)+(l-1)*wm));...
            yd;...
            zd-((zs+(hm/2)+(k-1)*hm))];
    end 
end

%To access R_kl's x,y,z coordinates, R_kl{k,l}(1/2/3)

R_kl_S = cellfun(@minus,R_kl,S,'Un',0);
R_kl_S_norm = cellfun(@(x)norm(x,2), R_kl_S);

R_kl_D = cellfun(@minus,R_kl,D,'Un',0);
R_kl_D_norm = cellfun(@(x)norm(x,2), R_kl_D);


for l = 1:nm
    for k = 1:nm
    R_kl_S_hat{k,l} = R_kl_S{k,l}./R_kl_S_norm(k,l);
    R_kl_D_hat{k,l} = R_kl_D{k,l}./R_kl_D_norm(k,l);
    end
end


for l = 1:nm
    for k = 1:nm
        N_kl_hat{k,l} = (R_kl_S_hat{k,l}+R_kl_D_hat{k,l})...
            /(2+2*transpose(R_kl_S_hat{k,l})*R_kl_D_hat{k,l})^1/2;
    end
end


for l = 1:nm
    for k = 1:nm
     beta{k,l} = asin(transpose(N_kl_hat{k,l})*[0;0;1]);
     alpha{k,l} = asin(...
         transpose(N_kl_hat{k,l})*[1;0;0]/cos(beta{k,l}));

    N_kl_hat{k,l} = [sin(beta{k,l})*cos(alpha{k,l});...
        cos(beta{k,l})*cos(alpha{k,l});...
        sin(alpha{k,l})];

    end
end


D_R_kl = cellfun(@minus,D,R_kl,'Un',0);
D_R_kl_norm = cellfun(@(x)norm(x,2), D_R_kl);

S_R_kl = cellfun(@minus,S,R_kl,'Un',0);
S_R_kl_norm = cellfun(@(x)norm(x,2), S_R_kl);


e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1];
e1_T = e1'; e2_T = e2'; e3_T = e3';
for l = 1:nm
    for k = 1:nm
        cos_theta_dash{k,l} = transpose(N_kl_hat{k,l})*(D_R_kl{k,l})...
            /D_R_kl_norm(k,l);
        R_kl_I = 2*cos_theta_dash{k,l}*N_kl_hat{k,l}-(D_R_kl{k,l})...
            /D_R_kl_norm(k,l);
        
        I_vec{k,l} = [e1_T*(R_kl{k,l}+e3_T*(S_R_kl{k,l})./(e3_T*R_kl_I)*R_kl_I);...
            e2_T*(R_kl{k,l}+e3_T*(S_R_kl{k,l})./(e3_T*R_kl_I)*R_kl_I);...
            e3_T*S{k,l}];
    end
end

for l = 1:nm
    for k = 1:nm
        theta_D_R{k,l} = e3_T*(R_kl{k,l}-I_vec{k,l})...
            /R_kl_D_norm(k,l);
    end
end

phi_half = 60;
m = -log(2)/log(cosd(phi_half));

I_mat = [1 0 0; 0 1 0; 0 0 1];

etta = 0.44;
resp = 0.54;
psi = deg2rad(70);
ref_indx = 1.5;
area = 10^-4;
for l = 1:nm
    for k = 1:nm
        I1 = (e1_T*S{k,l}-wm/2<=e1_T*I_vec{k,l}) * ...
            (e1_T*I_vec{k,l}<=e1_T*S{k,l}+wm/2);
        I2 = (e2_T*S{k,l}-hm/2<=e2_T*I_vec{k,l}) * ...
            (e2_T*I_vec{k,l}<=e2_T*S{k,l}+hm/2);
        E_kl{k,l} = (m+1)*rho/(2*pi)*cos_theta_dash{k,l}^(m)* ...
            e3_T*D_R_kl{k,l}*N_kl_hat{k,l}'*D_R_kl{k,l}./(D_R_kl_norm(k,l))^4 ...
            .*(1);
        
        if theta_D_R{k,l}<psi
             gain_PD_eve(k,l)=(ref_indx^2)/((sin(psi))^2);
        else
             gain_PD_eve(k,l)=0;   
        end
        
        h_IRS(k,l) = etta*resp*area*E_kl{k,l}*gain_PD_eve(k,l)*cos(theta_D_R{k,l});
            
%         h_los = 0; 
        A = 140e-3;
        
        source = [-0.26,2.5,0.5];
        dest = [xd(ind),yd,zd];
        h0 = 3;
        cos_m_11 = (h0./sqrt(norm(source-dest)).^(m+1));
        fov = deg2rad(70);
        n = 1.5;
        h_los = area*(m+1)*n^2*cos_m_11...
        ./(2*pi*(norm(source-dest))*(sin(fov))^2);
    
    
        R_sec(k,l) = 1/2*log((6*A^2*(h_IRS(k,l)+h_los)' *(h_IRS(k,l)+h_los)+3*pi*exp(1)*10^-18)/...
        (pi*exp(1)*A^2*(h_IRS(k,l)+h_los)' *(h_IRS(k,l)+h_los)+3*pi*exp(1)*10^-18));
    
    end
end

R_secret(ind) = max(max(R_sec));
h_IRS(k,l) = sum(sum(h_IRS));

end

h_0 = [6.08165852105811e-06,7.25772532383990e-06,8.59034867452604e-06,1.00589336060162e-05,1.16211253959985e-05,1.32098115518901e-05,1.47342191511320e-05,1.60867465443287e-05,1.71558069647752e-05,1.78428856412468e-05,1.80800080513214e-05,1.78428856412468e-05,1.71558069647752e-05,1.60867465443287e-05,1.47342191511320e-05,1.32098115518901e-05,1.16211253959985e-05,1.00589336060162e-05,8.59034867452604e-06,7.25772532383990e-06,6.08165852105811e-06];
h_opt = [6.97673770724660e-06,8.28325930015559e-06,9.73315792403789e-06,1.12901972993273e-05,1.28938495816673e-05,1.44587943391972e-05,1.58797220118320e-05,1.70423816139020e-05,1.78399132031617e-05,1.81912822126917e-05,1.80570926442438e-05,1.74481033404299e-05,1.64236612477904e-05,1.50803706415738e-05,1.35343606737961e-05,1.19021995462493e-05,1.02852165640842e-05,8.76010217092728e-06,7.37621445864377e-06,6.15810807462532e-06,5.11097483109363e-06];
h_rand = [4.86410586340751e-06,4.91197904451254e-06,6.94257437102338e-06,5.14162336596248e-06,5.84658775554916e-06,8.17739842282301e-06,8.01500976992924e-06,6.15878107922966e-06,1.25146915920678e-05,1.35672202083886e-05,1.16476171123668e-05,8.79296470211439e-06,1.01947299032549e-05,8.72536547113572e-06,5.94824623982233e-06,8.53915760537474e-06,6.18047172711482e-06,6.96626425825307e-06,4.46543492352523e-06,6.26431133744382e-06,4.78445797697000e-06];
y = [0.52,0.47,0.438,0.3960,0.3580,0.3340,0.3100,0.3,0.287,0.278,0.27,0,0,0,0.3100,0.3380,0.3580,0.4380,0.4700,0.5200];
x = -1:0.1:1;


figure('DefaultAxesFontSize',20);
plot(xd, h_opt,'linewidth',2);
hold on;
% plot(xd, h_0);
plot(xd, h_rand,'linewidth',2);
legend('Optimized','Random')
ylabel('Channel Gain') 
xlabel('Distance (m)')

figure('DefaultAxesFontSize',20);
plot(x(1:20),smooth(smooth(smooth(y))),'-*','linewidth',2)
hold on;
ylabel('Secrecy Rate'), xlabel('Distance between center of AP and Eve (m)');

