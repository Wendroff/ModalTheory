function output = select_script(run)

% This file contains Matlab codes for computing the eigenvalue
% problem corresponding to the 2D Helmholtz, and the 2D and 1D LNS
% equations. These are parts of the tutorial:
% "Modal Stability Theory" by Matthew Juniper, Ardeshir Hanifi,
% and Vassilios Theofilis, published in Applied Mechanics Reviews,
% 66(2), 2014.
%
% (c) 2014 vt_productions

% The main programs are
%
% Helmholtz2D.m    : solve the 2D Helmholtz EVP
% EVP2D            : set up and solve the 2d BiGlobal EVP
% LNSE1D           : solve the EVP pertinent to the 1d LNSE
%
% To execute, replace the argument 'run' by the string obtained
% from the corresponding function name, sans the suffix.

switch run
  case 'Helmholtz2D'
    output = Helmholtz2D();
  case 'EVP2D'
    output = EVP2D();
  case 'LNSE1D'
    output = LNSE1D();
  otherwise
    warning('AMR:nodemo', 'No such demo');
end

end

function output = Helmholtz2D()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Helmholtz2D:		Solve the 2D Helmholtz EVP                            %
%                                                                         %
%                       (c) 2014 vt_productions                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Set number of discretization nodes in the x and y spatial directions
Nx=48;
Ny=48;

% Set grids and collocation derivative matrices in the resolved spatial
% directions

Lx = 1;
Ly = 1;
[xgl,d1y,d2y] = set_chebyshev_derivative_matrices(Ny);
[xgrid,d1x,d2x] = map_y( 0., Lx, xgl, d1y, d2y);
[ygrid,d1y,d2y] = map_y( 0., Ly, xgl, d1y, d2y);

% Form the EVP and impose boundary conditions
[ A ] = form_Helmholtz_2D( Nx, Ny, d2x, d2y);
[ A ] = impose_Helmholtz_BCs( Nx, Ny, A );

% Solve the EVP using the QZ algorithm (full spectrum computation)
[V,lambda]=eig(A);

% Extract the eigenvalues
omegar=real(diag(lambda));
omegai=imag(diag(lambda));

% Plot the results
[phi] = plot_Helmholtz2D(-diag(lambda)/pi^2,V,ygrid,xgrid,Lx,Ly);

output = {};

end

function output = EVP2D()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% EVP2D:      Set up and solve the 2d BiGlobal EVP                        %
%                                                                         %
%                       (c) 2014 vt_productions                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

% Set number of discretization nodes in the x and y spatial directions
Nx=4;
Ny=128;

% Set flow parameters
Re = 10000;
alpha = 1;
% JHT test case:  OMEGAR=0.23752649, OMEGAI=0.00373967

% Re=7500;
% alpha=1;
% CHQZ test case: OMEGAR=0.24989154, OMEGAI=0.00223497

Lx=2*pi/alpha;

% Set grids and collocation derivative matrices in the resolved spatial
% directions
[xgrid]=set_x_grid(Nx);
[d1x,d2x]=set_fourier_derivative_matrices(Nx);
[xgrid,d1x,d2x] = map_x( 0., Lx, xgrid, d1x, d2x);
[xgl,d1y,d2y] = set_chebyshev_derivative_matrices(Ny);

% Form the EVP and impose boundary conditions
[ A, B ] = form_EVP2d( Nx, Ny, Re, xgl, d1x, d2x, d1y, d2y);
[ A, B ] = impose_BCs_EVP2d( Nx, Ny, d1y, A, B );

% Solve the EVP using the QZ algorithm (full spectrum computation)
[V,lambda]=eig(A,B,'qz');

% Extract the eigenvalues
omegar=real(diag(lambda));
omegai=imag(diag(lambda));

% Plot the results
plot_EVP2d(1i*diag(lambda),V,xgl,xgrid,Re,alpha)

output = {};

end

function output = LNSE1D()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% LNSE1D:      Solve the EVP pertinent to the 1d LNSE                     %
%                                                                         %
%                       (c) 2014 vt_productions                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% Set number of discretization nodes in the y spatial direction
Ny=128;

% Set flow and parameters
flow = 'SL';

% Set grid and collocation derivative matrices
[xgl,d1y,d2y] = set_chebyshev_derivative_matrices(Ny);

% If necessary, map the domain
[ y, d1y, d2y ] = mapping ( flow, xgl, d1y, d2y );

% Define the base flow
switch flow
  case 'PPF'
    
    [ U, Uy, Uyy ] = baseFlow_PPF( y );
    
    % Th94 test case:  OMEGAR=0.23752649, OMEGAI=0.00373967
    Re = 10000;
    alpha = 1;
    beta = 0;
    
    % CHQZ test case: OMEGAR=0.24989154, OMEGAI=0.00223497
    % Re=7500;
    % alpha=1;
    % beta=1;
    
  case 'SL'
    
    [ U, Uy, Uyy ] = baseFlow_SL( y );
    
    % Michalke case
    Re=1e5;
    alpha=0.4446;
    beta=0;
    
  otherwise
    
    disp('Flow Case Not Implemented')
    pause
end

% Form the EVP
[ A, B ] = form_LNSE1d(Re, alpha,  Ny, U, Uy, d1y, d2y);

% Impose BCs
[ A, B ] = impose_BCs_LNSE1d( Ny, A, B );

% Solve the EVP using the QZ algorithm (full spectrum computation)
[V,lambda]=eig(A,B,'qz');

% Extract the eigenvalues
omegar=real(diag(lambda));
omegai=imag(diag(lambda));

% Plot the results
plot_LNSE1d(flow,1i*lambda,V,y,Re,alpha,beta)

output = {};

end

function [xgl, d1y, d2y] = set_chebyshev_derivative_matrices(Ny)

% SET_CHEBYSHEV_DERIVATIVE_MATRICES: Define matrices for
%       Chebyshev spectral collocation differentiation
% INPUTS
% Ny:   number of Chebyshev Gauss-Lobatto collocation nodes
%
% OUTPUTS
% xgl:  Chebyshev Gauss-Lobatto points
% d1y:  1st CGL collocation derivative matrix
% d2y:  2nd CGL collocation derivative matrix
%
% (c) 2014 vt_productions

for j=0:Ny
  jj=j+1;
  xgl(jj)=cos(j*pi/Ny);
end

for j=0:Ny
  jj=j+1;
  if(j==0 | j==Ny )
    cj=2.;
  else
    cj=1.;
  end
  for k=0:Ny
    kk=k+1;
    if(k==0 | k==Ny )
      ck=2.;
    else
      ck=1.;
    end
    if(j~=k)
      d1y(jj,kk)=(cj/ck)*(-1)^(j+k)/(xgl(jj)-xgl(kk));
    end
  end
end
for j=1:Ny-1
  jj=j+1;
  d1y(jj,jj)=-xgl(jj)/(2.*(1.-xgl(jj)^2));
end
d1y(1,1)=(2.*Ny^2+1.)/6.;
d1y(Ny+1,Ny+1)=-d1y(1,1);

d2y = d1y*d1y;

end

function [ ygrid, d1y, d2y ] = map_y( yin, yout, ygrid, d1y, d2y )

% MAP_Y:           Linearly map a non-periodic domain
%
% INPUTS
% yin:             start of mapped domain
% yout:            end of mapped domain
%
% INPUTS/OUTPUTS
% ygrid, d1y, d2y: mapped domain (output)
%
% (c) 2014 vt_productions

a		= (yout - yin )/2.;
b		= (yin + yout )/2.;

ygrid   =  a * ygrid + b;

d1y     = 1./a * d1y;
d2y     = d1y * d1y;

end

function [ A ] = form_Helmholtz_2D( Nx, Ny, d2x, d2y)

% FORM_Helmholtz_2D: Form the 2D Laplacian
%
% INPUTS
% Nx, Ny:     number of Chebyshev collocation nodes
%
% d2x:        Chebyshev Gauss Lobatto collocation derivative matrix
% d2y:        Chebyshev Gauss Lobatto collocation derivative matrix
%
% OUTPUTS
% A:          LHS of the discretized Laplacian operator in 2D
%
% (c) 2014 vt_productions

Ix=eye(Nx+1);
Iy=eye(Ny+1);

A = kron(d2x,Iy)+kron(Ix,d2y);

end

function [ A ] = impose_Helmholtz_BCs( Nx, Ny, A )

% impose_Helmholtz_BCs: Impose homogeneous Dirichlet BCs
%
% INPUTS
% Nx, Ny:     number of Chebyshev collocation nodes
%
% INPUT/OUTPUT
% A:          LHS of the discretized Helmholtz operator in 2D
%
% (c) 2014 vt_productions

% .. phi=0 at y = Ly..

for i=1:Nx+1
  A(i*(Ny+1),:)=0;
  A(i*(Ny+1),i*(Ny+1))=1.;
end

% .. phi=0 at y = 0..

for i=1:Nx+1
  A(i+(i-1)*Ny,:)=0;
  A(i+(i-1)*Ny,i+(i-1)*Ny)=1.;
end

% .. phi=0 at x = Lx..

for i=2:Ny
  A(i,:)=0;
  A(i,i)=1.;
end

% .. phi=0 at x = 0..

for i=2:Ny
  A(Nx*(Ny+1)+i,:)=0;
  A(Nx*(Ny+1)+i,Nx*(Ny+1)+i)=1.;
end

end

function [phi]=plot_Helmholtz2D(eigval,eigvec,y,x,Lx,Ly)

% plot_Helmholtz2D plots the eigenfunctions computed by Helmholtz2D

Ny=length(y);
Nx=length(x);
NN=Nx*Ny;
[X,Y]=meshgrid(x,y);
x1=linspace(0,Lx,50);
y1=linspace(0,Ly,50);

[X1,Y1]=meshgrid(x1,y1);

disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
  
  figure(1)
  plot(real(eigval),imag(eigval),'o');
  axis([30,40,-1,1]);
  axis('square');
  ylabel('-imag(\lambda)/\pi^2')
  xlabel('-real(\lambda)/\pi^2')
  
  % pick the eigenfunction
  [xp,yp,button] = ginput(1);
  a=xp+sqrt(-1)*yp;
  [c,locate]=min(abs(eigval-a));
  
  % extract the eigenfunction
  phi=reshape(eigvec(1:NN,locate),Ny,Nx);
  
  % interpolate on a finer grid
  %    phi=interp2([X X],[Y Y],[phi phi],X1,Y1,'spline');
  
  phi=phi/max(max(phi));
  vecp=[0 0.2 0.4 0.6 0.8 1];
  figure(2)
  contour(X,Y,phi,vecp,'k-');
  axis('square');
  hold on;
  vecm=[-0.2 -0.4 -0.6 -0.8 -1];
  contour(X,Y,phi,vecm,'k-.');
  hold off;
end

end

function [ xgrid ] = set_x_grid( Nx )

% SET_X_GRID:   Define a periodic domain and discretize it
%               using Fourier spectral collocation
% INPUTS
% Nx:           number of Fourier collocation nodes
% OUTPUTS
% xgrid:        periodic domain
%
% (c) 2014 vt_productions

for jj=1:Nx
  j=jj-1;
  xgrid(jj)=2*pi*j/Nx;
end

end

function[d1x,d2x] = set_fourier_derivative_matrices(Nx)

% SET_FOURIER_DERIVATIVE_MATRICES: Define matrices for
%       Fourier spectral collocation differentiation
% INPUTS
% Nx:   number of Fourier collocation nodes
% OUTPUTS
% d1x:  1st Fourier collocation derivative matrix
% d2x:  2nd Fourier collocation derivative matrix
%
% (c) 2014 vt_productions

for jj = 1:Nx
  j=jj-1;
  for kk = 1:Nx
    k=kk-1;
    if(j==k)
      d1x(kk,jj)=0.;
      d2x(kk,jj)=-Nx^2/12-1./6;
    else
      d1x(kk,jj)=0.5*(-1)^(k+j)*cot((k-j)*pi/Nx);
      d2x(kk,jj)=-0.5*(-1)^(k+j)/sin((k-j)*pi/Nx)^2;
    end
  end
end

end

function [ xgrid, d1x, d2x ] = map_x( xin, xout, xgrid, d1x, d2x )

% MAP_X:           Linearly map a periodic domain
%
% INPUTS
% xin:             start of mapped domain
% xout:            end of mapped domain
%
% INPUTS/OUTPUTS
% xgrid, d1x, d2x: mapped domain (output)
%
% (c) 2014 vt_productions

L_x=(xout-xin);
metric  = 2.*pi/L_x;

xgrid   =  (L_x/2./pi) * xgrid + xin;

d1x     = metric * d1x;
d2x     = d1x * d1x;

end

function [ A, B ] = form_EVP2d( Nx, Ny, Re, xgl, d1x, d2x, d1y, d2y)

% FORM_EVP2D: Form the Linearized Navier-Stokes Equations in 2D
%
% INPUTS
% Nx:         number of Fourier collocation nodes
% d1x, d2x:   Fourier collocation derivative matrices
% Ny:         number of Chebyshev collocation nodes
% xgl:        Chebyshev Gauss Lobatto nodes
% d1y, d2y:   Chebyshev Gauss Lobatto collocation derivative matrices
% U, Uy:      Streamwise velocity of plane Poiseuille flow and 1st derivative
%
% OUTPUTS
% A:          LHS of the discretized LNSE operator in 2D
% B:          RHS of the discretized LNSE operator in 2D,
%             including an artificial compressibility parameter:
%             epsilon=1e-8
%
% (c) 2014 vt_productions

Ix=eye(Nx);
Iy=eye(Ny+1);

U=1.-xgl.^2;
Uy=-2*xgl;

epsilon = 1e-8;

Lop=1/Re*(kron(d2x,Iy)+kron(Ix,d2y))-kron(Ix,diag(U))*kron(d1x,Iy);

A11=Lop;
A12=-kron(Ix,diag(Uy));
A13=-kron(d1x,Iy);
A1=horzcat(A11,A12,A13);

A21=zeros(Nx*(Ny+1));
A22=Lop;
A23=-kron(Ix,d1y);
A2=horzcat(A21,A22,A23);

A31=kron(d1x,Iy);
A32=kron(Ix,d1y);
A33=zeros(Nx*(Ny+1));
A3=horzcat(A31,A32,A33);

A=vertcat(A1,A2,A3);

B11=eye(Nx*(Ny+1));
B12=zeros(Nx*(Ny+1));
B13=zeros(Nx*(Ny+1));
B1=horzcat(B11,B12,B13);

B21=zeros(Nx*(Ny+1));
B22=eye(Nx*(Ny+1));
B23=zeros(Nx*(Ny+1));
B2=horzcat(B21,B22,B23);

B31=zeros(Nx*(Ny+1));
B32=zeros(Nx*(Ny+1));
%B33=zeros(Nx*(Ny+1));
B33=epsilon*eye(Nx*(Ny+1));
B3=horzcat(B31,B32,B33);

B=vertcat(B1,B2,B3);

end

function [ A, B ] = impose_BCs_EVP2d( Nx, Ny, d1y, A, B )
%
% IMPOSE_BCS  Impose Boundary Conditions on the discretized LNSE operator.
% Impose BCs on perturbation velocity components only, none on pressure
%
% INPUTS/OUTPUTS
% A:          LHS of the discretized LNSE operator in 2D
% B:          RHS of the discretized LNSE operator in 2D
%
% (c) 2014 vt_productions

OffSet = Nx * (Ny+1);

% .. u=0 at y = +1..

for i=1:Nx
  A(i*(Ny+1)+0*OffSet,:)=0;
  B(i*(Ny+1)+0*OffSet,:)=0;
  A(i*(Ny+1)+0*OffSet,i*(Ny+1)+0*OffSet)=1.;
end

% .. u=0 at y = -1..

for i=1:Nx
  A(i+(i-1)*Ny+0*OffSet,:)=0;
  B(i+(i-1)*Ny+0*OffSet,:)=0;
  A(i+(i-1)*Ny+0*OffSet,i+(i-1)*Ny+0*OffSet)=1.;
end

% .. v=0 at y = +1..

for i=1:Nx
  A(i*(Ny+1)+1*OffSet,:)=0;
  B(i*(Ny+1)+1*OffSet,:)=0;
  A(i*(Ny+1)+1*OffSet,i*(Ny+1)+1*OffSet)=1.;
end

% .. v=0 at y = -1..

for i=1:Nx
  A(i+(i-1)*Ny+1*OffSet,:)=0;
  B(i+(i-1)*Ny+1*OffSet,:)=0;
  A(i+(i-1)*Ny+1*OffSet,i+(i-1)*Ny+1*OffSet)=1.;
end

end

function plot_EVP2d(eigval,eigvec,y,x,Re,alpha)

% plot_EVP2d plots the eigenfunctions computed by EVP2d

Ny=length(y);
Nx=length(x);
NN=Nx*Ny;
[X,Y]=meshgrid(x,y);
lambda=2*pi/alpha;
x1=linspace(0,x(end)+lambda,2*Nx*5);

disp('extend of x-domain'),x(end)
disp('# of points'),2*Nx*5

[X1,Y1]=meshgrid(x1,y);

disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
  
  figure(1)
  plot(real(eigval),imag(eigval),'o',[-100 100],[0 0],'r-');
  axis([0,1,-1,0.1]);axis('square');
  title(strcat({'Re= '},{num2str(Re)},{' , alpha= '},{num2str(alpha)}));
  ylabel('imag(\omega)')
  xlabel('real(\omega)')
  
  % pick the eigenfunction
  [xp,yp,button] = ginput(1);
  a=xp+sqrt(-1)*yp;
  [c,locate]=min(abs(eigval-a));
  
  % extract different components of eigenfunction
  u=reshape(eigvec(1:NN,locate),Ny,Nx);
  v=reshape(eigvec(NN+1:2*NN,locate),Ny,Nx);
  p=reshape(eigvec(2*NN+1:3*NN,locate),Ny,Nx);
  
  % interpolate on a finer grid
  u=interp2([X X+lambda],[Y Y],[u u],X1,Y1,'spline');
  v=interp2([X X+lambda],[Y Y],[v v],X1,Y1,'spline');
  p=interp2([X X+lambda],[Y Y],[p p],X1,Y1,'spline');
  
  % scale with respective maxima
  u = u/max(max(abs(u)));
  v = v/max(max(abs(v)));
  p = p/max(max(abs(p)));
  
  %define contours to plot
  vecp = [+0.1 +0.2 +0.3 +0.4 +0.5 +0.6 +0.7 +0.8 +0.9 +1];
  vecm = [-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1];
  
  % u
  figure(2)
  [C,h] = contour(X1/lambda,Y1,real(u),vecp,'k-');
  hold on;
  [C,h] = contour(X1/lambda,Y1,real(u),vecm,'k:');
  [C,h] = contour(X1/lambda,Y1,real(u),[0 0],'k.');
  hold off;
  axis('square');
  
  %    [C,h] = contourf(X1/lambda,Y1,real(u),10);
  %    set(h,'LineStyle','none')
  %    colormap('jet'); title('u'); ylabel('y'); xlabel('x/\lambda');
  %
  % v
  figure(3)
  [C,h] = contour(X1/lambda,Y1,real(v),vecp,'k-');
  hold on;
  [C,h] = contour(X1/lambda,Y1,real(v),vecm,'k:');
  [C,h] = contour(X1/lambda,Y1,real(v),[0 0],'k.');
  hold off;
  axis('square');
  
  %    [C,h] = contourf(X1/lambda,Y1,real(v),10);
  %    set(h,'LineStyle','none')
  %    colormap('jet'); title('v'); ylabel('y'); xlabel('x/\lambda');
  
  % Pressure
  figure(4)
  [C,h] = contour(X1/lambda,Y1,real(p),vecp,'k-');
  hold on;
  [C,h] = contour(X1/lambda,Y1,real(p),vecm,'k:');
  [C,h] = contour(X1/lambda,Y1,real(p),[0 0],'k.');
  hold off;
  axis('square');
  
  %    [C,h] = contourf(X1/lambda,Y1,real(p),10);
  %    set(h,'LineStyle','none')
  %    colormap('jet'); title('p'); ylabel('y'); xlabel('x/\lambda');
  
end

end

function [ y, d1, d2 ] = mapping ( flow_case, y, d1, d2 )

% mapping:         Map the [-1,1] CGL domain onto [-y_inf,y_inf]. 
%                  Compute the metric of the transformation and modify
%                  first and second derivative matrices accordingly.
%                  Plane Poiseuille flow and tanh(y) profiles implemented
%
% INPUT
% flow_case:       PPF and free shear layer (SL) model 
%
% INPUTS/OUTPUTS
% y, d1, d2:       CGL grid and collocation derivative matrices
%
% (c) 2014 vt_productions

metric(:,2) = zeros(size(y));

switch flow_case 
   
   case 'PPF'
        eta = y;
        metric(:,1) = ones(size(y));
        metric(:,2) = zeros(size(y));
   
   case 'SL'   
        len = 10;
        s_inf = 15;
        s = (len/s_inf)^2;
        for i=1:length(y)
            eta(i) = -len*y(i)/sqrt(1.+s-y(i)^2);
%       eta(i) = len*sinh(-s*y(i))/sinh(s);
        end
        for i=1:length(y)
            metric(i,1) = -sqrt(1.+s)*len^2/(len^2+eta(i)^2)^(3./2);
            metric(i,2) = 3*len^2*sqrt(1.+s)*eta(i)/(len^2+eta(i)^2)^(5./2);
%       metric(i,1) = -(1./len/s)*sinh(s)/sqrt(1.+(sinh(s)*eta(i)/len)^2);
%       metric(i,2) = (1./len^3/s)*eta(i)*sinh(s)^3/sqrt(1.+(sinh(s)*eta(i)/len)^2)^3;
        end
        
    otherwise
        disp('Flow Case Not Implemented')
        pause
end
     
y = eta;
for i=1:length(y)
    d1(i,:) = metric(i,1) * d1(i,:);
    d2(i,:) = metric(i,1)^2 * d2(i,:) + metric(i,2) * d1(i,:);
end
%d2 = d1*d1;
    
end

function [ U, Uy, Uyy ] = baseFlow_PPF( y )

% baseFlow_PPF: Define the parabolic velocity profile and its derivatives 
%               
% INPUTS
% y:            Chebyshev Gauss-Lobatto points

% OUTPUTS
% U, Uy, Uyy:   Streamwise velocity, first and second derivatives
%
% (c) 2014 vt_productions

U=1-y.^2;
Uy=-2*y;
Uyy=-2.;

end

function [ U, Uy, Uyy ] = baseFlow_SL( y )

% baseFlow_SL:  Define the hyperbolic tangent velocity profile and its derivatives 
%               
% INPUTS
% y:            Chebyshev Gauss-Lobatto points

% OUTPUTS
% U, Uy, Uyy:   Streamwise velocity, first and second derivatives
%
% (c) 2014 vt_productions

U=tanh(y);
for i=1:length(y)
    Uy(i)=1.-tanh(y(i))^2;
    Uyy(i)=-2.*tanh(y(i))*(1.-tanh(y(i))^2);
end

end

function [ A, B ] = form_LNSE1d(Re, alpha,  Ny, U, Uy, d1y, d2y)

% FORM_LNSE1D: Form the Linearized Navier-Stokes Equations in 1D

% INPUTS
% Re, alpha:  Reynolds number and streamwise wavenumber
% Ny:         number of Chebyshev collocation nodes
% d1y, d2y:   Chebyshev Gauss Lobatto collocation derivative matrices
% U, Uy:      Streamwise baseflow velocity and its 1st derivative 

% OUTPUTS
% A:          LHS of the discretized LNSE operator in 1D
% B:          RHS of the discretized LNSE operator in 1D,
%             including an artificial compressibility parameter: 
%             epsilon=1e-8

% (c) 2014 vt_productions

epsilon = 1e-8;

Lop=1/Re*d2y-alpha^2/Re*eye(Ny+1,Ny+1)-1i*alpha*diag(U);

A11=Lop;
A12=-diag(Uy);
A13=-1i*alpha*eye(Ny+1,Ny+1); 
A1=horzcat(A11,A12,A13);

A21=zeros(Ny+1);
A22=Lop;
A23=-d1y;
A2=horzcat(A21,A22,A23);

A31=1i*alpha*eye(Ny+1,Ny+1);
A32=d1y;
A33=zeros(Ny+1);
A3=horzcat(A31,A32,A33);

A=vertcat(A1,A2,A3);

B11=eye(Ny+1);
B12=zeros(Ny+1);
B13=zeros(Ny+1);
B1=horzcat(B11,B12,B13);

B21=zeros(Ny+1);
B22=eye(Ny+1);
B23=zeros(Ny+1);
B2=horzcat(B21,B22,B23);

B31=zeros(Ny+1);
B32=zeros(Ny+1);
%B33=epsilon*eye(Ny+1);
B33=zeros(Ny+1);
B3=horzcat(B31,B32,B33);

B=vertcat(B1,B2,B3);

end

function [ A, B ] = impose_BCs_LNSE1d( Ny, A, B )

% IMPOSE_BCS  Impose Boundary Conditions on the discretized LNSE operator.
% Impose BCs on perturbation velocity components only, none on pressure   

% INPUTS/OUTPUTS
% A:          LHS of the discretized LNSE operator in 1D
% B:          RHS of the discretized LNSE operator in 1D
%
% (c) 2014 vt_productions

OffSet = Ny+1;

% .. u=0 at y = +1..
A(Ny+1,:)=0;
B(Ny+1,:)=0;
A(Ny+1,Ny+1)=1.;

% .. u=0 at y = -1..
A(1,:)=0;
B(1,:)=0;
A(1,1)=1.;

% .. v=0 at y = +1..
A(Ny+1+1*OffSet,:)=0;
B(Ny+1+1*OffSet,:)=0;
A(Ny+1+1*OffSet,Ny+1+1*OffSet)=1.;

% .. v=0 at y = -1..
A(1+1*OffSet,:)=0;
B(1+1*OffSet,:)=0;
A(1+1*OffSet,1+1*OffSet)=1.;

end

function plot_LNSE1d(flow,eigval,eigvec,y,Re,alpha,beta)

%eglist=find(abs(eigval)<10)
%eigval=eigval(eglist);pause
%eigvec=eigvec(:,eglist)

N=length(y);

% Choose eigenvalue to plot

disp('Right click to pick the mode. Left click to stop');
button=1;
while button==1
    
    figure(1)
    plot(real(eigval),imag(eigval),'bo',[0 0.5],[0 0],'r-');
    axis('square');
    switch flow
        case 'PPF'
            ylim([-1,0.1]);
        case 'SL'
            ylim([-1,1]);
    end
    title(strcat({'Re= '},{num2str(Re)},{' , alpha= '},{num2str(alpha)},...
        {' , beta= '},{num2str(beta)}));
    ylabel('imag(\omega)')
    xlabel('real(\omega)')
    %ylim([-1,0.1]);
   
    % pick the eigenfunction    
    [xp,yp,button] = ginput(1);
    a=xp+sqrt(-1)*yp;
    [c,locate]=min(abs(eigval-a));

    % extract different components of eigenfunction
    u=eigvec(1:N,locate);
    v=eigvec(N+1:2*N,locate);
    p=eigvec(2*N+1:3*N,locate);
    
    u=u/max(max(abs(u)));
    v=v/max(max(abs(u)));
    p=p/max(max(abs(p)));
    
    figure(2)
    plot(abs(u),y,abs(v),y,abs(p),y);
    axis('square');
    legend('abs(u)','abs(v)','abs(p)')
    title(strcat({'Re= '},{num2str(Re)},{' , alpha= '},{num2str(alpha)},...
                 {', omega= '},{num2str(eigval(locate))}));
    ylabel('y')

end
dd = {y, abs(u), abs(v), abs(p)};
save eigfunc dd

end

function [ der1x, der2x ] = calculate_ddx( Nx, d1x, d2x, x )

%CALCULATE_DDX: Use Fourier spectral collocation to compute 
%               1st and 2nd derivatives of a periodic function
% INPUTS
% Nx:   number of (mapped) Fourier collocation nodes 
%       discretizing the periodic domain (not passed as an argument)
% x:    periodic function defined on the (mapped) domain
% d1x:  1st Fourier collocation derivative matrix, 
%       modified by the mapping transformation
% d2x:  2nd Fourier collocation derivative matrix, 
%       modified by the mapping transformation
% OUTPUTS
% der1x:1st derivative of the periodic function x
% der2x:2nd derivative of the periodic function x
%
% (c) 2014 vt_productions

for j =1:Nx
    der1x(j)=dot(d1x(j,:),x);
    der2x(j)=dot(d2x(j,:),x);
end
        
end

function [ der1y, der2y ] = calculate_ddy( Ny, d1y, d2y, y )

%CALCULATE_DDY: Use Chebyshev spectral collocation to compute 
%               1st and 2nd derivatives of a non-periodic function
% INPUTS
% Ny:   number of (mapped) Chebyshev collocation nodes 
%       discretizing the non-periodic domain (not passed as an argument)
% y:    periodic function defined on the (mapped) domain
% d1y:  1st Chebyshev collocation derivative matrix, 
%       modified by the mapping transformation
% d2y:  2nd Chebyshev collocation derivative matrix, 
%       modified by the mapping transformation
% OUTPUTS
% der1y:1st derivative of the non-periodic function x
% der2y:2nd derivative of the non-periodic function x


for j =1:Ny+1
    der1y(j)=dot(d1y(j,:),y);
    der2y(j)=dot(d2y(j,:),y);
end

end

function [ ygrid ] = set_y_grid( Ny )

% SET_Y_GRID:   Define the standard Chebyshev Gauss-Lobatto (CGL) 
%               spectral collocation points
% INPUTS
% Ny:           number of CGL collocation nodes
%
% OUTPUTS
% ygrid:        CGL domain

% (c) 2014 vt_productions

for jj=1:Ny+1
    j=jj-1;
    ygrid(jj)=cos(j*pi/Ny);
end

end