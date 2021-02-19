%CRJ330Config.m - Calculate the dimensional derivatives
%
% 22 Nov 2010
% INFORMATION - 
% Output: CRJ330Data:    Structure stores dimensional derivatives
%
% Original Author: Jason M.F. Zhang
% Revision: 1.0   $Date: 10/27/2009$


clear all

% ========== (1) flight conditions ==========
% (1.1) standard atmosphere &altitude
g = 32.2;              % ft/s^2, acceleration of gravity
alt = 33000;            % ft, altitude
rho = 7.9657e-004;      % slug/ft^3, air density

% (1.2) equilibrium
u0 = 726.5928;		    % ft/s, speed
theta0 = 0.0122;     % rad/s pitch angle
alpha0=0.0122;             % rad, AoA
gamma0 = 0;             % rad, flight path angle = 0
m = 1.2422e3;	        % slug, plane mass
J =[ 55717           0      -17789;
           0      369830           0;
      -17789           0      411017]; % slug*ft^2 martix of inertial
Ix=J(1,1);
Iy=J(2,2);
Iz=J(3,3);
Ixz=J(1,3);
Izx=J(3,1);
% ========== (2) Aircraft Geometric Configuration ==========
% (2.1) Wing
Wing.b = 67.8469;        % Wingspan (ft)
Wing.Cbar = 8.3208;	     % Wing mean geometric (aerodynamic) cord (ft)
Wind.Area= 520;          % Wing area (ft^2)
S=Wind.Area;
c_bar=Wing.Cbar;
b=Wing.b;
%Cw0=0; [TODO] - check which Cw shall be used
Cw0 = m*g / (0.5 * rho * u0^2 * S);
% ========== (3) Derivatives ===============================
% (3.1) Longitudinal
        CLo = 0.3654;
        CDo = 0.0352;
        Cmo = 0.0526;
        CTo = 0.0352;
        
       CmT= -0.0526;
       CLu= 0.2768;
       CDu= 0.0692;
       Cmu= 0.0335;
       CTu= 0;
      CmTu= 0;
       CLa= 6.8449;
       CDa= 0.3112;
       Cma= -4.1140;
      CmTa= 0;
    CLadot= 4.0521;
   
    CDadot= 0;
    Cmadot= -14.8779;
       CLq= 16.2541;
       CDq= 0;
       Cmq= -32.4358;
      CLde= 0.5926;
      CDde= 0;
      Cmde= -2.8814;
      
% Lateral
       Cyb= -1.0653;
       Clb= -0.0776;
       Cnb= 0.3123;
       Cyp= -0.2267;
       Clp= -0.6719;
       Cnp= -0.0016;
       Cyr= 0.6232;
       Clr= 0.2015;
       Cnr= -0.2788;
       Cyda= 0;
       Clda= 0.2210;
       Cnda= -0.0195;
       Cydr= 0.3099;
       Cldr= 0.0249;
       Cndr= -0.1392;
       
% Nondimensional Stability Derivatives in Body Axis
	
	% Longitudinal
	Cxu     = -(CDu + 2 * CDo) + ( CTu + 2 *  CTo);
	Cxa     = ( CLo -  CDa);
	Cxq     = - CDq;
	Cxa_dot = - CDadot;
	Czu     = -( CLu + 2 *  CLo);
	Cza     = -( CLa +  CDo);
	Czq     = - CLq;
	Cza_dot = - CLadot;
	Cmu     = ( Cmu  + 2 *  Cmo) + ( CmTu + 2 *  CmT);
	Cma     =  Cma;
	Cmq     =  Cmq;
	Cma_dot =  Cmadot;
	Cxde    = - CDde;
	Cxdp    =  CDo;
	Czde    = - CLde;
	Czdp    = 0;
	Cmde    =  Cmde;
	Cmdp    = 0;
	
	% Lateral
	Cyb     =  Cyb;
	Cyp     =  Cyp;
	Cyr     =  Cyr;
	Clb     =  Clb;
	Clp     =  Clp;
	Clr     =  Clr;
	Cnb     =  Cnb;
	Cnp     =  Cnp;
	Cnr     =  Cnr;
	Cyda    =  Cyda;
	Cydr    =  Cydr;
	Clda    =  Clda;
	Cldr    =  Cldr;
	Cnda    =  Cnda;
	Cndr    =  Cndr;

    
% Dimensional Stability Derivatives, Etkin's Formulation (Tables 4.4 and 4.5)
	
	% Longitudinal
	CRJData.DimStableDer.Longitudinal.Xu     = rho*u0*S*Cw0*sin(theta0) + 0.5*rho*u0*S*Cxu;     % <slug/s>
	CRJData.DimStableDer.Longitudinal.Xw     = 0.5*rho*u0*S*Cxa;                                % <slug/s>
	CRJData.DimStableDer.Longitudinal.Xq     = 0.25*rho*u0*c_bar*S*Cxq;                         % <slug*ft/s>
    CRJData.DimStableDer.Longitudinal.Xu_dot = 0;
    CRJData.DimStableDer.Longitudinal.Xw_dot = 0.25*rho*c_bar*S*Cxa_dot;                        % <slug>
    CRJData.DimStableDer.Longitudinal.Xq_dot = 0;
	CRJData.DimStableDer.Longitudinal.Zu     = -rho*u0*S*Cw0*cos(theta0) + 0.5*rho*u0*S*Czu;    % <slug/s>
	CRJData.DimStableDer.Longitudinal.Zw     = 0.5*rho*u0*S*Cza;                                % <slug/s>
	CRJData.DimStableDer.Longitudinal.Zq     = 0.25*rho*u0*c_bar*S*Czq;                         % <slug*ft/s>
    CRJData.DimStableDer.Longitudinal.Zu_dot = 0;
	CRJData.DimStableDer.Longitudinal.Zw_dot = 0.25*rho*c_bar*S*Cza_dot;                        % <slug>
    CRJData.DimStableDer.Longitudinal.Zq_dot = 0;
	CRJData.DimStableDer.Longitudinal.Mu     = 0.5*rho*u0*c_bar*S*Cmu;                          % <slug*ft/s>
	CRJData.DimStableDer.Longitudinal.Mw     = 0.5*rho*u0*c_bar*S*Cma;                          % <slug*ft/s>
	CRJData.DimStableDer.Longitudinal.Mq     = 0.25*rho*u0*c_bar^2*S*Cmq;                       % <slug*ft^2/s>
    CRJData.DimStableDer.Longitudinal.Mu_dot = 0;
	CRJData.DimStableDer.Longitudinal.Mw_dot = 0.25*rho*c_bar^2*S*Cma_dot;                      % <slug*ft>
    CRJData.DimStableDer.Longitudinal.Mq_dot = 0;
	CRJData.DimStableDer.Longitudinal.Xde    = 0.5*rho*u0^2*S*Cxde;                             % <slug*ft/s^2>
	CRJData.DimStableDer.Longitudinal.Xdp    = 0.5*rho*u0^2*S*Cxdp;                             % <slug*ft/s^2>
	CRJData.DimStableDer.Longitudinal.Zde    = 0.5*rho*u0^2*S*Czde;                             % <slug*ft/s^2>
	CRJData.DimStableDer.Longitudinal.Zdp    = 0.5*rho*u0^2*S*Czdp;                             % <slug*ft/s^2>
	CRJData.DimStableDer.Longitudinal.Mde    = 0.5*rho*u0^2*S*c_bar*Cmde;                       % <slug*ft^2/s^2>
	CRJData.DimStableDer.Longitudinal.Mdp    = 0.5*rho*u0^2*S*c_bar*Cmdp;                       % <slug*ft^2/s^2>

	% Lateral
	CRJData.DimStableDer.Lateral.Yv     = 0.5*rho*u0*S*Cyb;                                % <slug/s>
    CRJData.DimStableDer.Lateral.Yp     = 0.25*rho*u0*b*S*Cyp;                             % <slug*ft/s>  
    CRJData.DimStableDer.Lateral.Yr     = 0.25*rho*u0*b*S*Cyr;                             % <slug*ft/s> 
    CRJData.DimStableDer.Lateral.Yv_dot = 0;
    CRJData.DimStableDer.Lateral.Yp_dot = 0;
    CRJData.DimStableDer.Lateral.Yr_dot = 0;
    CRJData.DimStableDer.Lateral.Lv     = 0.5*rho*u0*b*S*Clb;                              % <slug*ft/s>
    CRJData.DimStableDer.Lateral.Lp     = 0.25*rho*u0*b^2*S*Clp;                           % <slug*ft^2/s>
    CRJData.DimStableDer.Lateral.Lr     = 0.25*rho*u0*b^2*S*Clr;                           % <slug*ft^2/s>
    CRJData.DimStableDer.Lateral.Lv_dot = 0;
    CRJData.DimStableDer.Lateral.Lp_dot = 0;
    CRJData.DimStableDer.Lateral.Lr_dot = 0;
    CRJData.DimStableDer.Lateral.Nv     = 0.5*rho*u0*b*S*Cnb;                              % <slug*ft/s>
    CRJData.DimStableDer.Lateral.Np     = 0.25*rho*u0*b^2*S*Cnp;                           % <slug*ft^2/s>
    CRJData.DimStableDer.Lateral.Nr     = 0.25*rho*u0*b^2*S*Cnr;                           % <slug*ft^2/s>
    CRJData.DimStableDer.Lateral.Nv_dot = 0;
    CRJData.DimStableDer.Lateral.Np_dot = 0;
    CRJData.DimStableDer.Lateral.Nr_dot = 0;
    CRJData.DimStableDer.Lateral.Yda    = 0.5*rho*u0^2*S*Cyda;                             % <slug*ft/s^2>
    CRJData.DimStableDer.Lateral.Ydr    = 0.5*rho*u0^2*S*Cydr;                             % <slug*ft/s^2>
    CRJData.DimStableDer.Lateral.Lda    = 0.5*rho*u0^2*S*b*Clda;                           % <slug*ft^2/s^2>
    CRJData.DimStableDer.Lateral.Ldr    = 0.5*rho*u0^2*S*b*Cldr;                           % <slug*ft^2/s^2>
    CRJData.DimStableDer.Lateral.Nda    = 0.5*rho*u0^2*S*b*Cnda;                           % <slug*ft^2/s^2>
    CRJData.DimStableDer.Lateral.Ndr    = 0.5*rho*u0^2*S*b*Cndr;                           % <slug*ft^2/s^2>
    
    
 CRJData.Config.m=m;
 CRJData.Config.g=g;
 CRJData.Config.J=J;

% CRJData.m=m;
% CRJData.g=g;
% CRJData.J=J;

    
    
    