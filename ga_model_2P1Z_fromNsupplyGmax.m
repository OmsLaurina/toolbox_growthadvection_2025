function varargout = ga_model_2P1Z_fromNsupplyGmax(Nsupply, Gmax1, Gmax2, varargin)


%% GA_MODEL_2P1Z_FROMNSUPPLY: plankton model
% Default parameters are based on Messié & Chavez (2017).
% 
% Use:
% [output=]ga_\_big from GA run initialized on May 1st, 2008model_2P2Z_fromNsupply(Nsupply,varargin)
%
% output is a structure containing 
%	.time, .Nsupply, .P_small, .P_big, .Z_small, .Z_big, 
%	.Nnew, .Nreg, .Chl, .PP, .u_small, .u_big, .g_small, .g_big, 
%	.units, .attributs
%
% Required input:
% 	Nsupply expressed in mmolC/m3/d, all units are carbon-based. Note - Nnew and Nreg represent new and regenerated nutrients, respectively
%		(NO3 and NH4 as a simplification) but are termed "Nnew" and "Nreg" to limit confusion with the unit, since they are expressed in carbon.
% 	Nsupply can be either a number, corresponding to the rate observed during upw_duration (default 1 day following Messié & Chavez 2017)
%					or a vector (then time needs to be given), that provides Nsupply as a function of time along a current trajectory,
%					for instance (useful to take Ekman pumping into account)
%
% Optional inputs:
% 'nbdays_advec'	number of days during which the model is run
% 'dt'				time step
% 'time'			can replace nbdays_advec and dt (required if Nsupply is a vector)
% 'upw_duration'	number of days during which Nsupply happens (default 1 day, not used if Nsupply is a vector)
% 'plot'			displays the plankton model outputs as a function of time
%
% Examples of use:
% ga_model_2P2Z_fromNsupply(1.3/16*106,'plot')				to reproduce Fig. 2 in Messié & Chavez 2017 (Nsupply = 1.3 mmolN/m3/d)
% ga_model_2P2Z_fromNsupply(11.2,'gmax_big',0.6*0.6,'eZ',0.1*0.6,'mZ',0.05*16/106*0.6,'plot')	to reproduce part of Fig. 1 in Messié et al. (in prep)
% 
% Monique Messié, 2021 for public version
% Reference: Messié, M., & Chavez, F. P. (2017). Nutrient supply, surface currents, and plankton dynamics predict zooplankton hotspots 
%					in coastal upwelling systems. Geophysical Research Letters, 44(17), 8979-8986, https://doi.org/10.1002/2017GL074322
% Differences with Messié and Chavez (2017):
%		all Zsmall excretion is now availabe as regenerated nutrients (ie no export on Zsmall excretion)
%		Zbig grazing formulation is different with the half-saturation constant applying to Z_small+P_big in both cases


%% -------------- Default parameters (see Messié & Chavez, 2017 suppl inf)

default_parameters = {...
        'umax1', 1.9872,   'umax2', 2.7648,  'kP1', 1,          'kP2', 3, ...
        'kZ1', 5,          'kZ2', 20,         'mP1', 0.1, ...
        'mP2', 0.2,        'm1Z', 0.1,        'm2Z', 0.061, ...
        'gamma', 0.6,      'epsilonP', 1,     'epsilon1Z', 0.3, ...
        'epsilon2Z', 0.7,  'P1_ini', 0.6,     'P2_ini', 0.1, ...
        'Z_ini', 0.6,      'PO4_ini', 0.5 ...
    };

[arg]=ga_read_varargin(varargin,[{'nb_days_advec',11,'dt',0.1,'time',[]},default_parameters]);
if length(Nsupply)>1 && isempty(arg.time), error('Give time if Nsupply is a vector'), end


%% -------------- Time

if isempty(arg.time), time=(time0:arg.dt:arg.nb_days_advec)'; 
else, time=arg.time(:); arg.dt=time(2)-time(1); 
end
nb_time=length(time); 

%% -------------- Nsupply

if length(Nsupply)==1
	Nsupply_max=Nsupply; 
	Nsupply=zeros(nb_time,1);
	Nsupply(:)=Nsupply_max; 
end

if length(Gmax1)==1
	Gmax1_max=Gmax1; 
	Gmax1=zeros(nb_time,1);
	Gmax1(:)=Gmax1_max; 
end

if length(Gmax2)==1
	Gmax2_max=Gmax2; 
	Gmax2=zeros(nb_time,1);
	Gmax2(:)=Gmax2_max; 
end

%% -------------- Initial conditions

P1 = time*NaN; P2 = time*NaN; Z = time*NaN; PO4 = time*NaN; u1 = time*NaN; u2 = time*NaN; g1 = time*NaN; g2 = time*NaN;
    PP1 = time*NaN; PP2 = time*NaN; G1 = time*NaN; G2 = time*NaN;
    exc = time*NaN; d_P1 = time*NaN; d_P2 = time*NaN; d1_Z = time*NaN; d2_Z = time*NaN;
    rec_P1 = time*NaN; rec_P2 = time*NaN; rec_Z = time*NaN; rec_exc = time*NaN;
    w_P1 = time*NaN; w_P2 = time*NaN; w1_Z = time*NaN; w2_Z = time*NaN; Export = time*NaN;


P1(1) = arg.P1_ini;
P2(1) = arg.P2_ini;
Z(1) = arg.Z_ini;
PO4(1) = arg.PO4_ini;



%% -------------- Loop on time

for t=2:nb_time

        % Taux de croissance (fonction de Monod)
        u1(t) = PO4(t-1) / (arg.kP1 + PO4(t-1)) * arg.umax1;
        u2(t) = PO4(t-1) / (arg.kP2 + PO4(t-1)) * arg.umax2;
        
        % Taux de consommation (fonction de Holling II)
        g1(t) = P1(t-1) / (arg.kZ1 + P1(t-1) + P2(t-1)) * Gmax1(t-1);
        g2(t) = P2(t-1) / (arg.kZ2 + P1(t-1) + P2(t-1)) * Gmax2(t-1);
        
        % Flux
        PP1(t) = u1(t) * P1(t-1);
        PP2(t) = u2(t) * P2(t-1);
        
        % GRAZING
        G1(t) = g1(t) * Z(t-1);
        G2(t) = g2(t) * Z(t-1);
        
        % Excrétion
        exc(t) = (1 - arg.gamma) * g1(t) * Z(t-1) + (1 - arg.gamma) * g2(t) * Z(t-1);
        
        % Mortalité
        d_P1(t) = arg.mP1 * P1(t-1);
        d_P2(t) = arg.mP2 * P2(t-1);
        d1_Z(t) = arg.m1Z * Z(t-1);
        d2_Z(t) = arg.m2Z * Z(t-1)^2;
        
        % Recyclage
        rec_P1(t) = arg.epsilonP * d_P1(t);
        rec_P2(t) = arg.epsilonP * d_P2(t);
        rec_Z(t) = arg.epsilon1Z * d1_Z(t);
        rec_exc(t) = arg.epsilon2Z * exc(t);
        
        % Lessivage
        w_P1(t) = (1 - arg.epsilonP) * d_P1(t);
        w_P2(t) = (1 - arg.epsilonP) * d_P2(t);
        w1_Z(t) = (1 - arg.epsilon1Z) * d1_Z(t);
        w2_Z(t) = (1 - arg.epsilon2Z) * exc(t);
        
        % Exportation
        Export(t) = w_P1(t) + w_P2(t) + w1_Z(t) + w2_Z(t) + d2_Z(t);
        
        % Mise à jour des concentrations de PO4 et biomasses
        max_available_PO4 = PO4(t-1) + Nsupply(t) * arg.dt + rec_P1(t) * arg.dt + rec_P2(t) * arg.dt + rec_Z(t) * arg.dt + rec_exc(t) * arg.dt;
        
        if PP1(t) > max_available_PO4 / arg.dt
            PP1(t) = max_available_PO4 / arg.dt;
        end
        
        max_available_PO4 = max_available_PO4 - PP1(t) * arg.dt;
        
        if PP2(t) > max_available_PO4 / arg.dt
            PP2(t) = max_available_PO4 / arg.dt;
        end
        
        % Biomasses
        P1(t) = P1(t-1) + PP1(t) * arg.dt - G1(t) * arg.dt - d_P1(t) * arg.dt;
        P2(t) = P2(t-1) + PP2(t) * arg.dt - G2(t) * arg.dt - d_P2(t) * arg.dt;
        Z(t) = Z(t-1) + arg.gamma * G1(t) * arg.dt + arg.gamma * G2(t) * arg.dt - d1_Z(t) * arg.dt - d2_Z(t) * arg.dt;
        PO4(t) = PO4(t-1) + Nsupply(t) * arg.dt + rec_P1(t) * arg.dt + rec_P2(t) * arg.dt + rec_Z(t) * arg.dt + rec_exc(t) * arg.dt - PP1(t) * arg.dt - PP2(t) * arg.dt;
        
end

%% -------------- Ouputs

units=struct('time','days','Nsupply','mmolC m^{-3} d^{-1}', 'Gmax1','d^{-1}', 'Gmax2','d^{-1}',...
	'P1','mmolC m^{-3}','P2','mmolC m^{-3}','Z','mmolC m^{-3}',...
	'PO4','mmolC m^{-3}', 'u1','d^{-1}','u2','d^{-1}','g1','d^{-1}','g2','d^{-1}');
output=struct('units',units,'time',time,'Nsupply',Nsupply, 'Gmax1', Gmax1, 'Gmax2', Gmax2,...
	'P1',P1,'P2',P2,'Z',Z,'PO4',PO4,...
	'u1',u1,'u2',u2,'g1',g1,'g2',g2,'attributs',struct('arg',arg));
varargout={output}; varargout=varargout(1:nargout);

return