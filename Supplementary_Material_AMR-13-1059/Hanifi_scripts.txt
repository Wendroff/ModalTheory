function output = select_script(run)

% This file contains Matlab codes for computing the eigenvalue
% problem corresponding to the Orr-Sommerfeld and linearized
% Navier-Stokes equations. These are parts of the tutorial:
% "Modal Stability Theory" by Matthew Juniper, Ardeshir Hanifi,
% and Vassilios Theofilis, published in Applied Mechanics Reviews,
% 66(2), 2014.
%
% (c) Ardeshir Hanifi 2014

% The main programs are
%
% run_OS.m    : run the temporal 1D stability code
% run_LNS.m   : run the temporal 1D stability problem
% run_PSE.m   : run the Parabolized Stability Equations
%
% Save the file as 'select_script.m'.
% To execute replace the argument 'run' by the string obtained
% from the corresponding function name, sans the suffix.

switch run
  case 'run_OS'
    output = run_OS();
  case 'run_LNS'
    output = run_LNS();
  case 'run_PSE'
    output = run_PSE();
  otherwise
    warning('AMR:nodemo', 'No such demo');
end

end

function output = run_OS()

% This script runs the temporal 1D stability code
%
% (c) Ardeshir Hanifi 2014

clear all
close all
%% Numerical parameters
N = 200; % number of Cehebychev modes

%% Flow parameters
Re=1e4;
alpha=1;
beta=0.4;

%% Select the collocation nodes and set up derivative operator
M=4; % Highest order of derivative required
[y, D] = chebdif(N, M);

%% Velocity profile corresponding to Plane Poiseuille Flow
U.u=1-y.^2;
U.uy=-2*y;
U.uyy=-2+y*0;
%% Set up the eigenvalue problem
[A,B] = OS_temp(U, D, alpha, beta, Re);

%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B);
eigval=diag(eigval);

%% Plot spectrum and eigenvectors
header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
  {', beta= '},{num2str(beta)});
plot_OS(eigval,eigvec,y,header)

output = {};

end

function output = run_LNS()

% This script runs the temporal 1D stability problem
%
% (c) Ardeshir Hanifi 2014

clear all
close all
%% Numerical parameters
N = 200; % number of Cehebychev modes

%% Flow parameters
Re=1e4;
alpha=1;
beta=0.4;

%% Select the collocation nodes and set up derivative operator
[y, D] = chebdif(N, 2);

%% Velocity profile corresponding to Plane Poiseuille Flow
U.u=1-y.^2;
U.uy=-2*y;
U.w=zeros(size(y)); % This field is required as the code is written for a 3D flow
U.wy=zeros(size(y)); % This field is required as the code is written for a 3D flow

%% Set up the eigenvalue problem
[A,B] = LNS_temp(U, D, alpha, beta, Re);

%% Solve the eigenvalue problem
[eigvec, eigval] = eig(A,B);
eigval=diag(eigval);

%% Plot spectrum and eigenvectors
header=strcat({'Re= '},{num2str(Re)},{', alpha= '},{num2str(alpha)},...
  {', beta= '},{num2str(beta)});
plot_LNS(eigval,eigvec,y,header)

output = {};

end

%------------------------------------------------------------------
function []=plot_OS(eigval,eigvec,y,header)
% []=plot_OS(eigval,eigvec,y,header)
% This function plots the spectrum and the corresponding eigenvectors
%
% INPUT
%   eigval: array containing the eigenvalues
%   eigvec: matrix corresponding the eigenvectors
%   y: normal coordinates
%   header: string to appear as the header of plot
%
% (c) Ardeshir Hanifi 2014

eglist=find(abs(eigval)<10);
eigval=eigval(eglist);
eigvec=eigvec(:,eglist);

% Choose eigenvalue to plot
disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
  
  figure(1)
  plot(real(eigval),imag(eigval),'o',[0 1],[0 0],'r-');
  ylim([-1,0.1]);
  title(header);
  ylabel('imag(\omega)')
  xlabel('real(\omega)')
  
  [xp,yp,button] = ginput(1);
  a=xp+sqrt(-1)*yp;
  [c,locate]=min(abs(eigval-a));
  
  v=eigvec(:,locate);
  figure(2)
  plot(real(v),y,imag(v),y);
  legend('real(v)','imag(v)')
  title(strcat(header,{', omega= '},{num2str(eigval(locate))}));
  ylabel('y')
  
end

end

%------------------------------------------------------------------
function []=plot_LNS(eigval,eigvec,y,header)
% []=plot_LNS(eigval,eigvec,y,header)
% This function plots the spectrum and selected eigenfunctions
%
% INPUT
%   eigval: array containing the eigenvalues
%   eigvec: matrix corresponding the eigenvectors
%   y: normal coordinates
%   header: string to appear as the header of plot
%
% (c) Ardeshir Hanifi 2014
%
eglist=find(abs(eigval)<10);
eigval=eigval(eglist);
eigvec=eigvec(:,eglist);

N=length(y);

% Choose eigenvalue to plot

disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
  
  figure(1)
  plot(real(eigval),imag(eigval),'o',[0 1],[0 0],'r-');
  ylim([-1,0.1]);
  title(header);
  ylabel('imag(\omega)')
  xlabel('real(\omega)')
  
  ylim([-1,0.1]);
  [xp,yp,button] = ginput(1);
  a=xp+sqrt(-1)*yp;
  [c,locate]=min(abs(eigval-a));
  
  u=eigvec(1:N,locate);
  v=eigvec(N+1:2*N,locate);
  w=eigvec(2*N+1:3*N,locate);
  p=eigvec(3*N+1:4*N,locate);
  
  figure(2)
  plot(abs(u),y,abs(v),y,abs(w),y,abs(p),y);
  legend('abs(u)','abs(v)','abs(w)','abs(p)')
  title(strcat(header,{', omega= '},{num2str(eigval(locate))}));
  ylabel('y')
  
end

end

%------------------------------------------------------------------
function [A,B] = OS_temp(U, D, alpha, beta, Re)
% [A,B]=OS_temp(U,D,ALPHA,BETA,RE)
% OS_temp sets up the operator corrsponding to local temporal 1D stability
% (Orr-Sommerfeldt) equation.
%
% INPUT
%   U: structure array containing meanflow quantities U.u, U.uy, U.uyy
%   D: spectral ifferential operator
%   alpha: stremwise wavenumber
%   beta: spanwise wavenumber
%   Re: Reynolds number
%
% OUTPUT
%   A and B: matrices corresponding to the eigenvalue problem A.q = omega.B.q
%
% (c) Ardeshir Hanifi 2014

N=length(U.u);
D1=D(:,:,1);
D2=D(:,:,2);
D4=D(:,:,4);
%------------------------------------------------------------------
% Set up operator matrices

zi = sqrt(-1);
k2 = alpha*alpha + beta*beta;
I = eye(N);

A=(D4-2*k2*D2)*zi/Re+diag(alpha*U.u)*(D2-k2*I)-diag(alpha*U.uyy);
B=D2-k2*I;
%------------------------------------------------------------------
% Boundary conditions
% Replace the empty rows of B with a multiple of same rows in A. This
% rmoves the singularity of the system and introduce eigenvalues which are
% located far from the physical ones by approperiate choice of multiple,
eps = 1e-4*zi;
Ov = zeros(1,N);

% v(ymax)=0
A(1,:) = Ov;
A(1,1) = 1;
B(1,:) = A(1,:)*eps;

% dv/dy(ymax)=0
A(2,:) = D1(1,:);
B(2,:) = A(2,:)*eps;

% dv/dy(0)=0
A(N-1,:) = D1(N,:);
B(N-1,:) = A(N-1,:)*eps;

% v(0) = 0
A(N,:) = Ov;
A(N,N) = 1;
B(N,:) = A(N,:)*eps;

end

%------------------------------------------------------------------
function [A,B] = LNS_temp(U, D, alpha, beta, Re)
% LNS_temp sets up the operator of general eigenvalue problem A.q = omega.B.q
% corrsponding to local temporal 1D stability equations.
%
% [A,B]=LNS_tmp(U,D,ALPHA,BETA,RE).
% 
% Input:
%   U:  structure containig streamwise and spanwise velocities and their 
%       first normal derivatives U.u, U.uy, U.w, U,wy
%   D:  sepctral differential operator
%   ALPHA: streamwise wavenumber
%   BETA: spanwise wavenumber
%   RE: Reynolds number
%
% OUTPUT
%   A and B: matrices corresponding to the eigenvalue problem A.q = omega.B.q
%
% (c) Ardeshir Hanifi, 2014

N=length(U.u);
D1=D(:,:,1);
D2=D(:,:,2);
%------------------------------------------------------------------
% Set up operator matrices

zi = sqrt(-1);
I = eye(N);
O = zeros(N);
xi=(D2-(beta^2+alpha^2)*I)/Re - zi*diag(alpha*U.u + beta*U.w);

A11 = zi*alpha*I;
A12 = D1;
A13 = zi*beta*I;
A21 = xi;
A22 = -diag(U.uy);
A24 = -zi*alpha*I;
A32 = xi;
A34 = -D1;
A42 = -diag(U.wy);
A43 = xi;
A44 = -zi*beta*I;

A = [A11, A12, A13, O  ; ...
     A21, A22, O,   A24; ...
     O,   A32, O,   A34; ...
     O,   A42, A43, A44];

%------------------------------------------------------------------

B = [O,     O,   O,     O; ...
     -zi*I, O,   O,     O; ...
     O,   -zi*I, O,     O; ...
     O,     O,   -zi*I, O];

%------------------------------------------------------------------
% Boundary conditions

Ov = zeros(1,4*N);

% u
A(N+1,:) = Ov;
A(N+1,1) = 1;
B(N+1,:) = Ov;

A(2*N,:) = Ov;
A(2*N, N) = 1;
B(2*N,:) = Ov;

%v
A(2*N+1,:) = Ov;
A(2*N+1, N+1) = 1;
B(2*N+1,:) = Ov;

A(3*N,:) = Ov;
A(3*N, 2*N) = 1;
B(3*N,:) = Ov;

%w
A(3*N+1,:) = Ov;
A(3*N+1, 2*N+1) = 1;
B(3*N+1,:) = Ov;

A(4*N,:) = Ov;
A(4*N, 3*N) = 1;
B(4*N,:) = Ov;

end

%------------------------------------------------------------------
function [x, DM] = chebdif(N, M)
%  The function [x, DM] =  chebdif(N,M) computes the differentiation 
%  matrices D1, D2, ..., DM on Chebyshev nodes. 
% 
%  Input:
%  N:        Size of differentiation matrix.        
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced 
%  accuracy suggested by W. Don and S. Solomonoff in 
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric 
%  identities to avoid the computation of differences 
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:   "Spectral Differencing with a Twist", by
%  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 

%  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by 
%  JACW, May 2003.

     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

     n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

     k = [0:N-1]';                        % Compute theta vector.
     th = k*pi/(N-1);

     x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

     T = repmat(th/2,1,N);                
     DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity. 
     DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
     DX(L) = ones(N,1);                    % Put 1's on the main diagonal of DX.

     C = toeplitz((-1).^k);               % C is the matrix with 
     C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
     C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
     Z(L) = zeros(N,1);                      % with zeros on the diagonal.

     D = eye(N);                          % D contains diff. matrices.
                                          
     for ell = 1:M
          D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
       D(L) = -sum(D');                            % Correct main diagonal of D
       DM(:,:,ell) = D;                                   % Store current D in DM
     end
end

%------------------------------------------------------------------
function output = run_PSE()
% pse.m is the main program for solving the PSE
% pse.m is the main program to solbe the Parabolized Stability Equations (PSE)
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%
%clear all
%% --------------------------------------------------------------------------
% Numerical Parameters
N=100;              % Number of nodes in wall-normal direction
ymax=70;           % Height of last node in wall-normal direction
stab=false;        % if stab=true the stablization of PSE is active
dpdx=1;            % if this parameter is =0 the x-derivative of pressure disturbance is neglected
FDorder='first';   % order of backward Euler discretization scheme (can be 'first' or 'second') 
iaux=1;            % auxiliary condition is based integral norm if (iaux=1) else on u_max
%% --------------------------------------------------------------------------
% Flow parameters
Re=400;            % Reynolds number at first x-station
Re_end=1000;       % Reynolds number at final x-station
Nstation=80;       % Number of x-stations
mfac=0.0;            % Falkner-Skan-Cook patrameter (U=x^mfac)
sweep=0;           % Sweep angle (in degree) at first x-station

%% --------------------------------------------------------------------------
% Perturbation parameters

beta=0.0;     % Spanwise wavenumber
F=1e-4;     % Reduced frequency
omega=F*Re; % Angular frequency

%% --------------------------------------------------------------------------
% Spectral discretisation operators (Chebycheff)

[y,D1,D2,W] = chebmat(N,ymax);

%% --------------------------------------------------------------------------
% Generate meanflow

disp('Computing mean flow')
x = linspace(Re,(Re_end/Re)^(2/(mfac+1))*Re,Nstation);
U = MeanFlow(x,y,Re,mfac,sweep);

%% --------------------------------------------------------------------------
% Solve local Orr-Sommerfeld eigenvalue problem to use as initial solution

alpha0=omega/0.36; % start approximation for TS waves
[q0,alpha0] = localOSS(alpha0,U(1),D1,D2,beta,omega,Re);

% Normalise initial disturbance to energy 1
q0 = q0/(0.5*sqrt(q0(1:3*N)'*blkdiag(W,W,W)*q0(1:3*N)));

% Plot the eigenfunction
c=['k','r','b','g'];
figure(1);clf(1);hold on; for k=1:4;plot(abs(q0((k-1)*N+1:k*N)),y,c(k));end
legend('u','v','w','p')

%% --------------------------------------------------------------------------
% Integrate the Parabolized Stability Equations

disp('Integrate Parabolized Stability Equations')

[q, alpha, ncalc] = integratePSE(q0, alpha0, x, y, U, D1, D2, W,...
    beta, omega, Re, Nstation, stab, dpdx, FDorder, iaux);
x=x(1:ncalc);

% Compute growth rate
sigma=growth( x,alpha,q,W );

% Compute N-factor
Nfactor=cumtrapz(x,sigma);


% Plot growth rate (scaled globally) and N-factor as function of local Re
Rex=(x/x(1)).^((mfac+1)/2)*Re;
figure(2);clf(2); plot(Rex,sigma);
grid on; xlabel('Re'); ylabel('growth rate');
figure(3);clf(3); plot(Rex,Nfactor);
grid on; xlabel('Re'); ylabel('N-factor');

d=([Rex' sigma]);
save('sigma.dat','d','-ascii','-double')

% %% --------------------------------------------------------------------------
% %  
% % This part of code is an example of how the local stability results could
% % be obtained.
% beta=0;     % Spanwise wavenumber
% F=1e-4;     % Reduced frequency
% 
% Rex=(400:20:800);
% alpha_l=zeros(size(Rex));
% for k=1:length(Rex)
% omega=F*Rex(k); % Angular frequency
% alpha0=omega/0.36; % start approximation for TS waves
% [q0,alpha_l(k)] = localOSS(alpha0,U(1),D1,D2,beta,omega,Rex(k));
% end
% 
% % plot results scaled by local refernce values
% figure(3);clf(3); plot(Rex,-imag(alpha_l));grid on;
%  
% %%
% % Data corresponding the neutral stability curve for Blasius flow
%   DNS_u=([340.4470  164.5750
%   329.7300  177.4480
%   323.2400  186.8020
%   317.8840  195.2120
%   313.7030  202.6430
%   310.5860  209.1890
%   308.1760  215.1470
%   306.0960  220.8290
%   304.1160  226.4280
%   302.3810  231.8230
%   301.0970  236.8430
%   300.4680  241.3180
%   300.5690  245.1850
%   301.2620  248.5580
%   302.4060  251.5560
%   304.1600  254.0470
%   306.9560  255.6700
%   311.2380  256.0560
%   317.4490  254.8360
%   326.0330  251.6400
%   337.4720  246.0680
%   352.6870  237.3510
%   371.0640  226.0020
%   388.6310  215.3270
%   403.5960  206.8190
%   416.9720  199.6340
%   429.6790  193.0060
%   441.8390  186.8330
%   453.4040  181.1550
%   464.3360  176.0050
%   474.6660  171.3550
%   484.4630  167.1500
%   493.7930  163.3330
%   502.7240  159.8480
%   511.3230  156.6400
%   519.6570  153.6520
%   527.7930  150.8290
%   535.8000  148.1140
%   543.7430  145.4510
%   551.6910  142.7850
%   559.7100  140.0600]);
% DNS_v=([344.1080  141.2510
%   308.5680  180.1410
%   303.3240  188.9460
%   294.3010  201.5050
%   288.2350  211.1260
%   284.4440  218.4880
%   281.1580  225.3490
%   277.9410  232.1410
%   274.9910  238.6690
%   272.4850  244.7550
%   270.4170  250.4060
%   268.6960  255.7130
%   267.2280  260.7680
%   265.9630  265.6220
%   264.9510  270.2250
%   264.2570  274.5110
%   263.9500  278.4140
%   264.0880  281.8750
%   264.6310  284.9330
%   265.4500  287.7170
%   266.4150  290.3570
%   267.5720  292.8050
%   269.4370  294.5510
%   272.2990  295.3070
%   276.0330  295.1960
%   280.4850  294.3730
%   285.7990  292.6940
%   295.1960  286.9600
%   309.6490  276.2050
%   325.5660  263.9960
%   340.5280  252.7360
%   354.0170  242.9390
%   366.6280  234.0130
%   378.5010  225.8200
%   389.4760  218.5200
%   399.4880  212.1750
%   408.6910  206.6330
%   417.2660  201.7150
%   425.3880  197.2480
%   433.1590  193.1290
%   440.6300  189.3070
%   447.8520  185.7330
%   454.8730  182.3590
%   461.7440  179.1330
%   468.4980  176.0240
%   475.1410  173.0250
%   481.6740  170.1350
%   488.0990  167.3520
%   494.4180  164.6740
%   500.6340  162.1000
%   506.7480  159.6260
%   512.7620  157.2510
%   518.6780  154.9740
%   524.4980  152.7920
%   530.2240  150.7040
%   535.8580  148.7060
%   541.4030  146.7980
%   546.8590  144.9780
%   552.2290  143.2430
%   557.5150  141.5910
%   562.7190  140.0210]);
% Exp_data=([430.9300  100.0000
%   711.4530  100.0000
%   337.3260  162.0000
%   473.8370  162.0000
%   310.9300  200.0000
%   405.0580  200.0000
%   284.4190  250.0000
%   324.8260  250.0000]);
% plot(DNS_u(:,1),DNS_u(:,2),DNS_v(:,1),DNS_v(:,2),Exp_data(:,1),Exp_data(:,2),'o')
% legend('DNS (u) (Berlin et al. 1998)','DNS (v) (Berlin et al. 1998)', 'Exp. (Klingmann et al 1993)')
% ylim([0,300]);xlabel('Re');ylabel('F');title('Neutral stability curve for Blasius flow');
% grid on

output = {};

end

%------------------------------------------------------------------
function auxval = auxfunc(eulco, y, q1, q2, q3, W, iaux)
% auxval = auxfunc(eulco, q1, q2, q, W)
% This function computes the value of the auxiliary function 
% coresponding to integral of conjg(q)*dq/dx.
%
%   eulco: contains the dsicretization coefficients
%   y: values of normal coordinates
%   q1, q2 and q3: the eigenfunctions at stations i-2, i-1 and i-th 
%                 (the currenet) streamwise station
%   W: contains spectral inetgral weights for Chebuychev polynomials
%   iaux: =1 integral norm, =2 max(u)
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%


N=length(y);

% Compute umax
umax1 = findqmax(y,q1(1:N));
umax2 = findqmax(y,q2(1:N));
umax3 = findqmax(y,q3(1:N));

% Compute the auxiliary condition
if (iaux==1)
    % Compute derivative of q
    dqdx = eulco(1)*q1 + eulco(2)*q2 + eulco(3)*q3;
    auxval = q3'*blkdiag(W,W,W,W)*dqdx;
else
    auxval = (eulco(1)*umax1 + eulco(2)*umax2 + eulco(3)*umax3);
end

end

%------------------------------------------------------------------
function [y,D1,D2,W] = chebmat(N,ymax)
% [z,D1,D2,W] = chebmat(N,zmax)
% This function sets up the integration weights and differential operatores 
% corresponding to the spectral discretization using Chenychev polynomials. 
% The differential operators are also transormed using a mapping from 
% the spectral domain [-1,1] to [0,ymax].
%
% N: order of Chebychem polynomial is N-1
% ymax: Height of last node in wall-normal direction
%
% The output variable are
% y: coordinates of Gauss-Lobbato nodes in physical space [0,ymax]
% D1, D2: first and second derivative operatires
% W: spectral integraation weights
%
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

[yc, D] = chebdif(N, 2);
D1 = D(:,:,1);
D2 = D(:,:,2);

% Define the mapping to the physical domain
y=(yc+1)*ymax/2;

% Map Derivative operators onto interval [0 ymax]
D1 = D1*(2/ymax);
D2 = D2*(2/ymax)^2;

% Get integration weights
W = iwt(N)*ymax/2;

end

%------------------------------------------------------------------
function [eulco] = compeulco(order, dx1, dx2)
%  [eulco] = compeulco(order, dx1, dx2) computes the coefficients for backward Euler FD scheme
%
% Input
%   order: 'first' or 'second' for first- or second-order backward Euler scheme
%   dx1: x(i)-x(i-1)
%   dx2: x(i-1)-x(i-2)
%
% Output:
%   eulco: an array containing the discretization coefficients such that 
%          df/dx=eulco(1)*f(i-2)+eulco(2)*f(i-1)+eulco(3)*f(i)
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

if strcmp(order, 'first')
    eulco=([0, -1, 1])/dx1;
elseif strcmp(order, 'second')
    d=dx2+dx1;
    eulco = ([dx1/(dx2*d), -d/(dx1*dx2), (2*dx1+dx2)/(dx1*d)]);
end

end
%------------------------------------------------------------------
function [ fx ] = dfx( f,x,difforder )
% [ fx ] = dfx( f,x,difforder )
% DFX Computes derivatives of F with respect to x
% 
% Input
%   f: funxtion values
%   x: values of the variable
%   difforder: =1 first-order Euler, =2 second-order Euler scheme
%
% Output
%   fx: df/dx
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

idim=length(f);
fx=zeros(size(f));
fx(1)=(f(1)-f(2))/(x(1)-x(2));
for i=2:idim
    if (min(difforder,i-1)==1)             % first order Euler
        dx1=x(i)-x(i-1);
        fx(i)=(f(i)-f(i-1))/dx1;
    else                          % second order Euler
        dx1=x(i)-x(i-1);
        dx2=x(i-1)-x(i-2);
        c2=dx1/(dx2*(dx2+dx1));
        c1=-(dx2+dx1)/(dx1*dx2);
        c0=(2.*dx1+dx2)/(dx1*(dx2+dx1));
        fx(i)=c2*f(i-2)+c1*f(i-1)+c0*f(i);
    end
end

end

%------------------------------------------------------------------
function [ x, x2, f2 ] = extrapol( x1, x2, f1, f2, iter)
% x = extrapol( x1, x2, f1, f2)
% This function find x value such that f(x)=0 using a secant method.
%
% Input:
%   f1 and f2: function values corresponding to variables x1 and x2
%   iter: index of iteration at curent iteration step
%
% Output:
%   x: estimated value of x correspopnding to f(x)=0
%   x2 and f2: same as in the input.
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

if (iter==1)
    x = x2*1.001;
else
    x = x2 - (x2 - x1)/(f2 - f1)*f2; 
end

end

%------------------------------------------------------------------
function [qmax] = findqmax(y,q)
% [qmax] = findqmax(q,y)
% This function finds the value of q with the maximum absolute value 
% by ploynomial fitting.
%
% Input:
%    y: array containing coordinates of grid points
%    q: function values at gridpoints given in y
%
% Output:
%    qmax: estimated value of q with maximum absolute value
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

N=length(y);

% Find the maximum of abs(q) on the original grid
qabs=abs(q);
[qmax, imax]=max(qabs);

if (imax==1 || imax==N)
    disp('Findqmax: max values at boundaries')
    qmax=q(imax);
    return;
end

% Fit polynomial to data
Np=2; % order of polynomial
i1=max(1,imax-floor(Np/2));
i2=min(i1+Np,N);
p = polyfit(y(i1:i2),qabs(i1:i2),Np);

% find zero of derivative of the polynomial
dp = polyder(p);
y0=[y(i1) y(i2)];
fun=@(x)polyval(dp,x);
ymax=fzero(fun,y0);

% Find the maximum on the interpolated grid
p = polyfit(y(i1:i2),q(i1:i2),Np);
qmax=polyval(p,ymax);

end

%------------------------------------------------------------------
function dy = fsc(t,y,mfac)
% dy = fsc(t,y,mfac) sets up the functions coresponding to
% Falkner-Skan-Cooke simililarity soloution
%
% Input:
%   t: variable required by ode45
%   y: array containing function values corresponding to differential equation for similarity solution
%   mfac: Falkne-Skan velocity parameter U=x^mfac 
% 
% Output:
%   dy: derivative of function y
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

dy=zeros(5,1);
betaH=2*mfac/(mfac+1);
dy(1)=y(2);
dy(2)=y(3);
dy(2)=y(3);
dy(3)=-(y(1)*y(3)+betaH*(1-y(2)*y(2)));
dy(4)=y(5);
dy(5)=-y(1)*y(5);
end

%------------------------------------------------------------------
function [ sigma ] = growth( x,alpha,q,W )
% [ sigma ] = growth( x,alpha,q,W )
% GROWTH computes the nonlocal growthrate of perturbations
%
% Input
%   x: streamwise coordinate of stations
%   alpha: vector containig streamwise wavenumber at all x-stations
%   q:  matrix containing eigenfunctions at all x-stations
%   W: is the weightfunction of spectral integration
%
% Output 
%    sigma: the nonlocal growth rate.
%
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

N=length(q(:,1))/4;

% Growthrate based on disturbance energy
energy=zeros(size(alpha));
umax=zeros(size(alpha));
for k=1:length(x)
    energy(k)=abs(0.5*q(1:3*N,k)'*blkdiag(W,W,W)*q(1:3*N,k));            
end
sigma=-imag(alpha)+dfx( log(sqrt(energy)),x,1 );
sigma(1)=-imag(alpha(1)); % Local growth rate

end

%------------------------------------------------------------------
function [q, alphax, ncalc] = integratePSE(q0, alpha0, x, y, U, D1, D2, W,...
                       beta, omega, Re, Nstation, stab, dpdx, FDorder, iaux)
% [q,alphax,ncalc] = integratePSE(q0, alpha0, x, y, U, D1, D2, W,...
%                           beta, omega, Re, Nstation, stab, dpdx, FDorder, iaux)                   
% Thuis function integrates the Parabolised Stability Equations (PSE).
%
% Input
%    q0, alpha0: complex eigenfunction and stermawise wavenumber at initial station
%    x, y: arrays containing streamwise and normal coordinates of the grid
%    U: structure array containing the meanflow quantities
%    D1, D2, W: first- and second derivative operators and spectral integration weights
%    beta, omega, Re: spanwise wavenumber, frequency and Reynolds number
%    Nstation: number of stations in x-direction to compute
%    stab: if stab=true the stablization of PSE is active
%    dpdx: if =0 the x-derivative of pressure disturbance is neglected
%    FDorder: order of backward Euler discretization scheme (can be 'first' or 'second') 
%    iaux: auxiliary condition is based integral norm if (iaux=1) else on u_max
%
% Output
%    q: complex array contatining eigen functions (u,v,w,p) at computed x-stations
%    alphax: complex array containing streamwise wavenumber at computed x-stations
%    ncalc: number of sucessfully computed stations
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

%initialise variables
N=length(y);
q=zeros(N*4,Nstation);        q(:,1) = q0;
alphax=zeros(Nstation,1);     alphax(1) = alpha0;
alpha = alpha0;
local=false;

dx= [0 diff(x)];
alpha1=0;
auxval1=0;

% integration loop
n=1;
ncalc=Nstation;
 while (n < Nstation) 
    n=n+1;
    disp(' ');
    disp(['Station: ', num2str(n)]);
    
    %decide which FDorder for backward scheme
    if (n == 2 || strcmp(FDorder, 'first'))
        scheme_order = 'first';
    elseif (n ~= 2 && strcmp(FDorder, 'second'))
        scheme_order = 'second';
    end

    %get coefficients for backward scheme
    [eulco] = compeulco(scheme_order, dx(n), dx(n-1));
    
    %loop over alpha to satisfy auxiliary condition
    iloop = 0;
    auxval = 1e8;

    while abs(auxval) > 1e-08 && (iloop<50)
        iloop = iloop + 1;

        %stablization coefficient
        s = 0;
        if stab
            if ( real(alpha) ~= 0  && dx(n) <= 1/abs(real(alpha))*1.01 )
                s = (0.5/abs(real(alpha)) - dx(n)*.5)*2;
            end
            %disp(['stabilisation: s = ' num2str(s)]);
        end
        
        %set up operator matrices with boundary conditions
        [RHS, LHS] = operator(local,U(n),alpha, D1, D2,...
                     beta,omega,Re,eulco,q(:,max(1,n-2)),q(:,n-1),s,dpdx);

        %solve equation system
        q(:,n) = LHS\RHS;
        
        %compute the auxiliary condition
        auxval = auxfunc(eulco, y, q(:,max(1,n-2)), q(:,n-1), q(:,n), W, iaux);
        disp(['It: ', num2str(iloop), ' ,alpha: ', num2str(alpha), ' , auxval: ', num2str(abs(auxval))]);
        alphax(n) = alpha;
        
        %estimate new alpha
        [ alpha, alpha1, auxval1 ] = extrapol( alpha1, alpha, auxval1, auxval, iloop);

    end
    if abs(auxval) > 1e-08
        ncalc=n-1;
        q=q(:,1:ncalc);
        alphax=alphax(1:ncalc);
        n=Nstation+100; 
    end
    
 end
end

%------------------------------------------------------------------
function W = iwt(N)
% IWT: sets up the weigth function for the spectral integration and
% form the matrix used for computation of disturbance energy.
%
% Input
%    N: order of Chebychev olynomial is N-1
%
% Output:
%    W: array containing the integration weights
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%


N1=N-1;
%compute chebycheff integration weights
n = 0:1:N1;
j = 0:1:N1;
b = ones(1, N);
b([1 N]) = 0.5;
c = 2*b;
b = b/N1;
S = cos(n(3:N)'*j*(pi/N1));
W = diag(b.*((2+(c(3:N).*((1+(-1).^n(3:N))./(1-n(3:N).^2)))*S)));

end

%------------------------------------------------------------------
function [q,alpha] = localOSS(alphainit, U, D1, D2, beta, omega, Re)
% [q,alpha] = localOSS(alphainit, U, D1, D2, beta, omega, Re) computes eigenvalue of 
% local spatial stability equation
%
% Input:
%    alphainit: start approximation of eigenvalue
%    U: structure containing mean flow velocity field and its derivatives
%    D1 and D2: matrices corresponding to operator of first- and second derivative
%    beta: spanwise wavenumber of perurbation
%    omega: angular frequency of perturbation
%    Re: Reynolds number
%
% Outout:
%   q: eigenfunction
%   alpha: eigenvalue (streamwise wavenumber)
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

% Find eigenvalue of spatial OSS system using a secant method
loop = true;
iloop = 0;
alpha = alphainit;
local=true;
N=length(U.u);

alpha1=0;
v01=0;

disp([' ']);
disp(['**** Local analysis for Re= ',num2str(Re),', F=',num2str(omega/Re)]);
disp([' ']);
while loop
    iloop = iloop +1;
    
    [RHS,LHS] = operator(local, U, alpha, D1, D2, beta, omega, Re);
    q = LHS\RHS;
    v02 = q(2*N);
    [ alpha, alpha1, v01 ] = extrapol( alpha1, alpha, v01, v02, iloop);
    
    disp(['It: ', num2str(iloop), ' ,alpha: ', num2str(alpha), ' , w0: ', num2str(abs(v02))]);
    
    if abs(v02) <= 2e-8
        loop = false;
    end
    
    if (iloop>50) error ('**** No converged solution found'); end

end

end

%------------------------------------------------------------------
function U = MeanFlow(x,y,Re,mfac,sweep)
% U = MeanFlow(x,z,Re,mfac,sweep)
% This function computes the mean flow quantities for a Falkner-Skan-Cooke 
% boundary layer. 
% 
% x: array containing coordinates of grid in x-direction
% y: array containing coordinates of grid in z-direction
% Re: Reynolds number at initial x-station
% mfac: Falkner-Skan-Cooke parameter U=C*x^mfac
% sweep: flowsweep angle at initial x-station
%
% Output is the structre array 'U' which contains 
% U(i).u, U(i).v, U(i).w: streamwise, normal and spanwise velocity
% U(i).uy, U(i).vy, U(i).wy: normal derivatives of u, v and w 
% U(i).ux, U(i).wx: streawwise derivatives of u and w 
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%

Nstations = length(x);

% generate self-semilar Falkner-Skan-Cooke profiles
Ufsc = similar_fsc(y,mfac);

% generate meanflow data at first station
we=tan(sweep*pi/180);
U(1).u = Ufsc.u;
U(1).v = Ufsc.v/Re;
U(1).w = Ufsc.w*we;
U(1).uy = Ufsc.uy;
U(1).wy = Ufsc.wy*tan(sweep*pi/180);
U(1).vy = Ufsc.vy/Re;
U(1).ux = Ufsc.uy.*( y*(mfac-1)*0.5 )/Re + mfac*Ufsc.u/Re;
U(1).wx = Ufsc.wy.*( y*(mfac-1)*0.5 )*we/Re;

% set global reference scales
method = 'spline';

% rescale meanflow such that all quantities are based on reference scales of
% first station and interpolate onto Gauss-Lobatto Grid
for i = 1:Nstations
    ue = (x(i)/x(1))^mfac;
    lscale = sqrt(x(i)/x(1)/ue);
    rescale = sqrt(ue*x(i)/x(1));
    yi=y*lscale;
    ve = U(1).v(1)*ue/rescale-U(1).ux(1)*ue/rescale*(yi(1)-y(1));
    
    U(i).u=interp1(yi,U(1).u*ue,y,method,ue);
    U(i).w=interp1(yi,U(1).w   ,y,method,we);
    U(i).v=interp1(yi,U(1).v*ue/rescale,y,method,ve);
 
    U(i).uy=interp1(yi,U(1).uy*ue/lscale,y,method,0);
    U(i).wy=interp1(yi,U(1).wy/lscale   ,y,method,0);
    U(i).vy=interp1(yi,U(1).vy*ue/lscale/rescale,y,method,U(1).vy(1)*ue/lscale/rescale);

    U(i).ux=interp1(yi,U(1).ux*ue/rescale^2,y,method,U(1).ux(1)*ue/rescale^2);
    U(i).wx=interp1(yi,U(1).wx/rescale^2   ,y,method,U(1).wx(1)/rescale^2);
    
end

end

%------------------------------------------------------------------
function [RHS,LHS] = operator(local, U, alpha, D1, D2, beta, omega, ...
                              Re, varargin)
% OPERATOR sets up the operator corrsponding to local and nonlocal stability
% equations.
%
% [RHS,LHS]=OPERATOR(LOCAL,U,ALPHA,D1,D2,BETA,OMEGA,RE,N,EULCO,Q1,Q2,S)
% If LOCAL=false parameters EULCO, Q1, Q2 and S are optional.
%
% Input
%	local: if =true the local operator corresponding to local stability theory is set up
%	U: structure containing meanflow field and its derivatives
%	alpha: streamwise wavenumber
%	D1 and D2: matrices corresponding to first- and second derivative operators
%	beta: spanwise wavenumber
%	omega: angular frequency
%	Re: Reynolds number
%	eulco: coefficients of Euler backward discretization scheme
%	q1 and q2: amplitudefunctions at x(i-1) and x(i-2)
%	s: coefficient of stabilization term
%   dpdx: if =0 the x-derivative of pressure disturbance is neglected
%
% Output
%	LHS and RHS: matrices corresponding to PSE, LHS.q=RHS 
%
% (c) Ardeshir Hanifi & David Tempelmann, 2014
%
%------------------------------------------------------------------
N=length(U.u);

if (local==true) 
    eulco = zeros(3,1);
    q1 = zeros(N*4,1);
    q2 = zeros(N*4,1);
    s = 0;
    dpdx = 0;
elseif (local==false) 
    if (nargin~=13)
        disp('**** OPERATOR: Wrong number of parameters')
        disp('[RHS,LHS] = operator(local,U,alpha,D1,D2,beta,omega,Re,eulco,q1,q2,s)')
        error('Stopping programe');
    else
        eulco = varargin{1};
        q1 = varargin{2};
        q2 = varargin{3};
        s = varargin{4};
        dpdx = varargin{5};
        if (dpdx~=0); dpdx=1; end
    end
end

%------------------------------------------------------------------
% Set up operator matrices
% A.q + B.q_y + C.q_yy + D.q_x=0

if (local==true); inoloc=0; else inoloc=1; end
i = sqrt(-1);
I = eye(size(D1));
O = zeros(size(D1));
xi=diag(-i*omega + (beta^2+alpha^2)/Re + i*(alpha*U.u + beta*U.w));

A11 = i*alpha*I;
A13 = i*beta*I;
A21 = xi+diag(U.ux)*inoloc;
A22 = diag(U.uy);
A24 = i*alpha*I;
A32 = xi+diag(U.vy)*inoloc;
A41 = diag(U.wx)*inoloc;
A42 = diag(U.wy);
A43 = xi;
A44 = i*beta*I;

A = [A11, O,   A13, O  ; ...
     A21, A22, O  , A24; ...
     O,   A32, O  , O  ;
     A41, A42, A43, A44]; ...

%------------------------------------------------------------------

VD1 = diag(U.v)*D1*inoloc;

B = [O,   D1,  O,   O;  ...
     VD1, O,   O,   O;  ...
     O,   VD1, O,   D1; ...
     O,   O,   VD1, O ];

%------------------------------------------------------------------

C = [O,      O,      O,      O; ...
     -D2/Re, O,      O,      O; ...
     O,      -D2/Re, O,      O; ...
     O,      O,      -D2/Re, O];

%------------------------------------------------------------------
if (local)   
    D=zeros(size(A));
else

    u = diag(U.u);
    
    D = [I, O, O, O     ; ...
         u, O, O, I*dpdx; ...
         O, u, O, O     ; ...
         O, O, u, O ];

    D = D + s*(A+B+C); %Include stabilisation term
end

%Left-hand-side and Right-hand-side of equation system
LHS = A + B + C + D*eulco(3);
RHS = -D*(eulco(1)*q1 + eulco(2)*q2);

%---------------------------------------------------------------------
%Boundary conditions
Ov = zeros(1,4*N);

% u(ymax)=0
LHS(N+1,:) = Ov;
LHS(N+1,1) = 1;
RHS(N+1) = 0;

% u(0)=0
LHS(2*N,:) = Ov;
LHS(2*N, N) = 1;
RHS(2*N) = 0;

% v(ymax)=0
LHS(2*N+1,:) = Ov;
LHS(2*N+1, N+1) = 1;
RHS(2*N+1) = 0;

if (local)
    % p(0)=1
    LHS(3*N,:) = Ov;
    LHS(3*N, 4*N) = 1;
    RHS(3*N) = 1.0;
else
    % v(0)=0
    LHS(3*N,:) = Ov;
    LHS(3*N, 2*N) = 1;
    RHS(3*N) = 0;
end

% w(ymax)=0
LHS(3*N+1,:) = Ov;
LHS(3*N+1, 2*N+1) = 1;
RHS(3*N+1) = 0;

% w(0)=0
LHS(4*N,:) = Ov;
LHS(4*N, 3*N) = 1;
RHS(4*N) = 0;

end

%------------------------------------------------------------------
function Ufsc=similar_fsc(y,mfac)
% Ufsc=SIMILAR_FSC(YIN,MFAC)
% This function computes the selfsimilar Falkner-Skan-Cook boundary-layer profiles
%
% Input
%	y:  contains normal coordinates weher velocity values are required.
%	mfac:  Falkner-Skan-Cooke parameter U=C*x^mfac. 
%
% Output 
%	Ufsc.u, Ufsc.v, Ufsc.w: streamwise, normal and spanwise velocity 
%	Ufsc.uy, Ufsc.vy, Ufsc.wy:  normal derivatives of u, v, and w
%
% (c) Ardeshir Hanifi & David Tempelmann
%
%-------------------------------------------------------------------------------
% Start approximation for f'' and g' 
fbis= -0.0791*mfac^4+0.7414*mfac^3-1.4302*mfac^2+1.6813*mfac+0.3318;
fbis=fbis*2/(mfac+1);
gprim= -1.3429*mfac^4+2.6745*mfac^3-1.8293*mfac^2+0.7102*mfac+0.3244;
gprim=gprim*2/(mfac+1);

% Options for ODE solver
options = odeset('RelTol',1e-12,'AbsTol',[1e-8 1e-8 1e-8 1e-8 1e-8]);
tspan=[0 max(y)];
ic=[0 0 fbis 0 gprim];

% Integrate equations until convergend solution is found
[yt,f]=ode45(@(t,y) fsc(t,y,mfac),tspan,ic,options);
ue=f(end,2);
ue1=ue;
fbis1=fbis;
fbis=fbis+0.01;
tol=1e-8;
while abs(ue-1)>tol
    ic=[0 0 fbis 0 1];
    [yt,f]=ode45(@(t,y) fsc(t,y,mfac),tspan,ic,options);
    ue=f(end,2);
    dfbis=(1-ue)*(fbis-fbis1)/(ue-ue1);
    fbis1=fbis;
    ue1=ue;
    fbis=fbis+dfbis;
end

% Interpolate the velocity profiles on the input grid
c=1/sqrt((mfac+1)/2);
vv=((1.0-mfac)*f(:,2).*yt-(1.0+mfac)*f(:,1))*.5*c;
vvz=((1.0-mfac)*f(:,2)+(1.0-mfac)*f(:,3).*yt-(1.0+mfac)*f(:,2))*.5*c;
Ufsc.u=interp1(yt*c,f(:,2),y,'cubic');
Ufsc.w=interp1(yt*c,f(:,4)/f(end,4),y,'cubic');
Ufsc.v=interp1(yt*c,vv,y,'cubic');
Ufsc.uy=interp1(yt*c,f(:,3)/c,y,'cubic');
Ufsc.wy=interp1(yt*c,f(:,5)/f(end,4)/c,y,'cubic');
Ufsc.vy=interp1(yt*c,vvz/c,y,'cubic');

end


