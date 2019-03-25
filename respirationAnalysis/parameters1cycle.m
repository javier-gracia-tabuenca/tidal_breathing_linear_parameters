function [ parValues , parNames, parValues_conf ] = parameters1cycle( Vin, Fin, Vex, Fex, fs, varargin )
% [ parValues , parNames, Vin, Fin, vrEx, frEx, h, parValues_conf ] = parameters1cycle( Vin, Fin, Vex, Fex, fs, varargin )
% Given a respiratory cycle spleted in inspiratory volume "Vin"; inspiratory flow Fin";
% expiratory volume "Vex"; and expiratory flow Fex"; with sampling frequency "fs" 
% calculates the linear parameters given in parNames
%
% Further information below before each algorithm/step
%
% INPUTS:
% Vin = inspiratory volume
% Fin = inspiratory flow
% Vex = expiratory volume
% Fex = expiratory flow
% fs = sampling frequency
%
% OUTPUTS:
% parValues = values of the calculated parametes in the order given in parNames
% parNames = label names of the caculated parameters
% parValues_conf = confident values of the parameter in the order given in
% parNames, if the calculation of a parameter does not involve a confidence
% value the defaut is 1
%
%
% VARARGIN
% 'normalizeV' = normalised volume to 1
% 'plot' = plot final result
% 'plot_color', c = set the color of the signal 
%
% PARAMETERS
% 'Vt' = tidal volume 
% 'tI' = inspiration time
% 'tE' = expiration time
% 'PTIF' = peak tidal inspiratory flow
% 't_PTIF' = time to peak tidal inspiratory flow
% 'MTIF' = mean tidal inspiratory flow
% 'tPTIF_per_tI' = time to peak tidal inspiratory flow / inspiration time 
% 'VPTIF_per_VI' = volume to peak tidal inspiratory flow / inpiration volume 
% 'TIF50' = flow at 50% inspired volume
% 'PTEF' = peak tidal expiratory flow 
% 't_PTEF' = time to peak tidal expiratory flow 
% 'MTEF' = mean tidal expiratory flow,
% 'tPTEF_per_tE' = time to peak tidal expiration flow / expiration time  
% 'V_PTEF' = 
% 'VPTEF_per_VE' = 
% 'TEF50' = flow at 50% remaining expiration volume
% 'TEF25' = flow at 25% remaining expiration volume
% previous parameters as defined in [1]
% 'S','Sintercept','S2' = [2]
% 'Trs','EV' = [3]
% 'SPvf' = [4]
%
%
% REFERENCES 
%
% [1]Beydon N, et al. 
% Thoracic Society/European Respiratory Society Statement: Pulmonary Function Testing
% in Preschool Children. Am J Respir Crit Care Med 175: 1304–1345, 2007.
% [2] Williams E, Madgwick R, Morris M. Tidal expired airflow patterns in
% adults with airway obstruction. Eur Respir J 12: 1118 –1123, 1998
% [3]Morris MJ, Madgwick RG, Collyer I, Denby F, Lane DJ. Analysis of expiratory tidal
% flow patterns as a diagnostic tool in airflow obstruction. Eur Respir J 12: 1113–1117,
% 1998.
% [4] Gracia-Tabuenca J, Tidal breathing flow profiles during sleep phases 
% in wheezing infants measured by impedance pneumography. 
% J Appl Physiol, 2019.
%
% Author: javier.gracia.tabuenca@gmail.com    Date: 12.12.2016
%

n = 0;
SD = 'none';
plotflag = '';
normalize = 0;
plot_ix = 0;
h_axes = [];
c='k';
c2='b';
while n < length(varargin)
	n = n + 1;
	if strcmp(varargin{n}, 'plot')
		plotflag ='plot';
	elseif strcmp(varargin{n}, 'parNames')
		n = n + 1;
		parNames=varargin{n};
	elseif strcmp(varargin{n}, 'plot_ix')
		n = n + 1;
		plot_ix=varargin{n};
	elseif strcmp(varargin{n}, 'plot_color')
		n = n + 1;
		c2=varargin{n};
		plotflag ='plot';
	elseif strcmp(varargin{n}, 'h_axes')
		n = n + 1;
		h_axes=varargin{n};
		plotflag='plot';
	elseif strcmp(varargin{n}, 'normalizeV')%not tested
		normalize=1;
	end
	
end



parValues = zeros(1,length(parNames));
parValues_conf = zeros(1,length(parNames));

h=[];
%normalize volume to 1
if normalize
	% remove volume trend
	v=[ Vin ; Vex ];
	trend = linspace(v(1),v(end),length(v))';
	v=v-trend;
	Vin = v(1:length(Vin));
	Vex = v(length(Vin)+1:end);
	
	%inspiration
	k = range(Vin);
	Vin = Vin ./ k;
	Fin = Fin ./ k;
	%expiration
	k = range(Vex);
	Vex = Vex ./ k;
	Fex = Fex ./ k;
end


% ploting variables
fsize=8;
v=[ Vin ; Vex ];
f=[ Fin ; Fex ];
vr= max(v)-v;
fr=-f;
t=getT(f,fs)-length(Fin)/fs;
if strcmp(plotflag,'plot')
	if isempty(h_axes)
		h(1) = subplot(3,2,1);
		h(2) = subplot(3,2,2);
		h(3) = subplot(3,2,3);
		h(4) = subplot(3,2,4);
		h(5) = subplot(3,2,5);
		h(6) = subplot(3,2,6);
	else
		h=h_axes;
	end
	marginV=0.1*abs(max(v)-min(v));
	marginF=0.1*abs(max(f)-min(f));
	marginT=0.1*abs(t(end)-t(1));
	
	textV = min(v);
	textF = min(f);
	%vol
	axes(h(1))
	hold on
	plot(t,v,'Color',c2);
	xlim([t(1)-marginT t(end)+marginT])
	ylim([min(v)-marginV max(v)+marginV])
	grid on
	title('volume')
	xlabel('t [s]')
	ylabel('vol [l] or [k_{\Omega]}')
	%flow
	axes(h(2))
	hold on
	plot(t,f,'Color',c2);
	xlim([t(1)-marginT t(end)+marginT])
	ylim([min(f)-marginV max(f)+marginV])
	grid on
	title('flow')
	xlabel('t [s]')
	ylabel('flow [l/s] or [k_{\Omega]}/s]')
	%F_V
	axes(h(3))
	plot(vr,fr,'Color',c2);
	hold on
	xlim([min(vr)-marginV max(vr)+marginV])
	ylim([min(fr)-marginV max(fr)+marginV])
	%FLIP FIGURE
	grid on
	title('FVloop')
	xlabel('vol [l] or [k_{\Omega]}')
	ylabel('flow [l/s] or [k_{\Omega]/s}]')
	%flow
	if length(h)>3
		axes(h(4))
		hold on
		grid on
		title('Normalize Pos-spiratory-peak flow')
		xlabel('t [%]' )
		ylabel('flow [%]')
	end
end

vrEx = vr(length(Vin)+1:end);
frEx = fr(length(Vin)+1:end);



plot_offset_vol=0.05*abs(max(v)-min(v));
plot_offset_flow=0.05*abs(max(f)-min(f));
%Vt
ix = find(strcmp(parNames,'Vt'));
if ix
	Vt = Vin(end)-Vin(1);
	parValues(ix) = Vt;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(1));
		line([t(length(Vin)) t(length(Vin))],[Vin(1) Vin(end)],'LineStyle','--','Color',c)
		plot(t(length(Vin)),Vin(end),[c 'o'])
		text(t(length(Vin)),Vin(end)+plot_ix*plot_offset_vol,sprintf('Vf=%.2f',Vt),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end


%IBI
ix = find(strcmp(parNames,'IBI'));
if ix
	T = length(v)/fs;
	parValues(ix) = T;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(1));
		text(t(round(length(v)/2)),min(v)-plot_ix*plot_offset_vol,sprintf('IBI=%.2fs',T),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end

%tI
ix = find(strcmp(parNames,'tI'));
if ix
	tI = length(Vin)/fs;
	parValues(ix) = tI;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(1));
		text(t(round(length(Vin)/2)),min(v)+plot_ix*plot_offset_vol,sprintf('tI=%.2fs',tI),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%tE
ix = find(strcmp(parNames,'tE'));
if ix
	tE = length(Vex)/fs;
	parValues(ix) = tE;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(1));
		text(t(round(length(Vin)+length(Vex)/2)),min(v)+plot_ix*plot_offset_vol,sprintf('tE=%.2fs',tE),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%PTIF
ix = find(strcmp(parNames,'PTIF'));
if ix
	[ptif,Mix]= max(Fin);
	parValues(ix) = ptif;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		line([t(Mix) t(Mix)],[0 ptif],'LineStyle','--','Color',c)
		plot(t(Mix),ptif,[c 'o'])
		text(t(Mix),ptif+plot_ix*plot_offset_flow,sprintf('PTIF=%.2f',ptif),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
		axes(h(3));
		line([vr(Mix) vr(Mix)],[0 -ptif],'LineStyle','--','Color',c)
		plot(vr(Mix),-ptif,[c 'o'])
		text(vr(Mix),-(ptif+plot_ix*plot_offset_flow),sprintf('PTIF=%.2f',ptif),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%PTEF
ix = find(strcmp(parNames,'PTEF'));
if ix
	[ptef,mix]= min(Fex);
	parValues(ix) = -ptef;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		line([t(length(Fin)+mix) t(length(Fin)+mix)],[0 ptef],'LineStyle','--','Color',c)
		plot(t(length(Fin)+mix),ptef,[c 'o'])
		text(t(length(Fin)+mix),ptef+plot_ix*plot_offset_flow,sprintf('PTEF=%.2f',-ptef),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		axes(h(3));
		line([vr(length(Fin)+mix) vr(length(Fin)+mix)],[0 -ptef],'LineStyle','--','Color',c)
		plot(vr(length(Fin)+mix),-ptef,[c 'o'])
		text(vr(length(Fin)+mix),-ptef+plot_ix*plot_offset_flow,sprintf('PTEF=%.2f',-ptef),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%TEF25_PTEF
ix = find(strcmp(parNames,'TEF25_PTEF'));
if ix
	[ptef,mix]= min(Fex);
	ptef = -ptef;
	[~,mix]=min(abs(Vex-0.25));
	tef25 = -Fex(mix);
	parValues(ix) = tef25/ptef;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		line([vr(length(Fin)+mix) vr(length(Fin)+mix)],[0 tef25],'LineStyle','--','Color',c)
		plot(vr(length(Fin)+mix),tef25,'o','Color',c)
		text(vr(length(Fin)+mix),tef25+plot_ix*plot_offset_flow,sprintf('TEF25/PTEF=%.2f',tef25),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%t_PTIF
ix = find(strcmp(parNames,'t_PTIF'));
if ix
	[ptif,Mix]= max(Fin);
	t_ptif=(length(Fin)-Mix)/fs;
	parValues(ix) = t_ptif;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		plot(t(Mix),0,[c 'o'])
		text(t(Mix),-plot_ix*plot_offset_flow,sprintf('t_{PTIF}=%.2f',t_ptif),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%t_PTEF
ix = find(strcmp(parNames,'t_PTEF'));
if ix
	[ptef,mix]= min(Fex);
	t_ptef=mix/fs;
	parValues(ix) = t_ptef;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		plot(t(length(Fin)+mix),0,[c 'o'])
		text(t(length(Fin)+mix),plot_ix*plot_offset_flow,sprintf('t_{PTEF}=%.2f',t_ptef),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%V_PTIF
ix = find(strcmp(parNames,'V_PTIF'));
if ix
	[ptif,Mix]= max(Fin);
	v_ptif = max(Vin)-Vin(Mix);
	parValues(ix) = v_ptif;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		plot(vr(Mix),0,'o','Color',c)
		text(vr(Mix),-plot_ix*plot_offset_flow,sprintf('V_{PTIF}=%.2f',v_ptif),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end
%V_PTEF
ix = find(strcmp(parNames,'V_PTEF'));
if ix
	[ptef,mix]= min(Fex);
	v_ptef = max(Vex)-Vex(mix);
	parValues(ix) = v_ptef;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		plot(vr(length(Fin)+mix),0,'o','Color',c)
		text(vr(length(Fin)+mix),-plot_ix*plot_offset_flow,sprintf('V_{PTEF}=%.2f',v_ptef),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end


%TIF50
ix = find(strcmp(parNames,'TIF50'));
if ix
	[~,mix]=min(abs(Vin-0.5));
	tif50 = Fin(mix);
	parValues(ix) = tif50;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		line([vr(mix) vr(mix)],[0 -tif50],'LineStyle','--','Color',c)
		plot(vr(mix),-tif50,'o','Color',c)
		text(vr(mix),-tif50-plot_ix*plot_offset_flow,sprintf('TIF50=%.2f',tif50),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%TEF50
ix = find(strcmp(parNames,'TEF50'));
if ix
	[~,mix]=min(abs(Vex-0.5));
	tef50 = -Fex(mix);
	parValues(ix) = tef50;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		line([vr(length(Fin)+mix) vr(length(Fin)+mix)],[0 tef50],'LineStyle','--','Color',c)
		plot(vr(length(Fin)+mix),tef50,'o','Color',c)
		text(vr(length(Fin)+mix),tef50+plot_ix*plot_offset_flow,sprintf('TEF50=%.2f',tef50),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end

%TEF25
ix = find(strcmp(parNames,'TEF25'));
if ix
	[~,mix]=min(abs(Vex-0.25));
	tef25 = -Fex(mix);
	parValues(ix) = tef25;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		line([vr(length(Fin)+mix) vr(length(Fin)+mix)],[0 tef25],'LineStyle','--','Color',c)
		plot(vr(length(Fin)+mix),tef25,'o','Color',c)
		text(vr(length(Fin)+mix),tef25+plot_ix*plot_offset_flow,sprintf('TEF25=%.2f',tef25),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',fsize)
	end
end


%MTIF
ix = find(strcmp(parNames,'MTIF'));
if ix
	mtif = mean(Fin);
	parValues(ix) = mtif;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		text(min(t),textF-plot_ix*plot_offset_flow,sprintf('MTIF=%.2f',mtif),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%MTEF
ix = find(strcmp(parNames,'MTEF'));
if ix
	mtef = mean(Fex);
	parValues(ix) = mtef;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		text(max(t),textF-plot_ix*plot_offset_flow,sprintf('MTEF=%.2f',mtef),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end


%tPTIF_per_tI
ix = find(strcmp(parNames,'tPTIF_per_tI'));
if ix
	[ptif,mix]= max(Fin);
	t_ptif=(length(Fin)-Mix)/fs;
	tI = length(Fin)/fs;
	parValues(ix) = t_ptif/tI;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		text(max(t)/4,textF-plot_ix*plot_offset_flow,sprintf('t_{PTIF}/t_I=%.2f',t_ptif/tI),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%tPTEF_per_tE
ix = find(strcmp(parNames,'tPTEF_per_tE'));
if ix
	[ptef,mix]= min(Fex);
	t_ptef=mix/fs;
	tE = length(Fex)/fs;
	parValues(ix) = t_ptef/tE;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(2));
		text(max(t)*3/4,textF-plot_ix*plot_offset_flow,sprintf('t_{PTEF}/t_E=%.2f',t_ptef/tE),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end



%VPTIF_per_VI
ix = find(strcmp(parNames,'VPTIF_per_VI'));
if ix
	[ptif,Mix]= max(Fin);
	v_ptif = max(Vin)-Vin(Mix);
	vi = Vin(end)-Vin(1);
	parValues(ix) = v_ptif/vi;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		text(min(v),textF-plot_ix*plot_offset_flow,sprintf('V_{PTIF}/V_I=%.2f',v_ptif/vi),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end
%VPTEF_per_VE
ix = find(strcmp(parNames,'VPTEF_per_VE'));
if ix
	[ptef,mix]= min(Fex);
	v_ptef = max(Vex)-Vex(mix);
	ve = Vex(1)-Vex(end);
	parValues(ix) = v_ptef/ve;
	parValues_conf(ix) = 1;
	if strcmp(plotflag,'plot')
		axes(h(3));
		text(max(v),textF-plot_ix*plot_offset_flow,sprintf('V_{PTEF}/V_E=%.2f',v_ptef/ve),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
	end
end



% PARAMETER S
if ~isempty(find(strcmp(parNames,'Sn'))) || ~isempty(find(strcmp(parNames,'Sintercept'))) || ~isempty(find(strcmp(parNames,'S2'))) || ~isempty(find(strcmp(parNames,'SE2')))
	[ptef,mix]= min(Fex);
	Fpp = -Fex(mix : end); % Post-peak section of flow signal
	if length(Fpp) == 1 %if there is a problewm i the signal
		Fpp = [0;Fpp;Fpp;10];
	end
	nFpp = Fpp-Fpp(end);
	nFpp = nFpp/nFpp(1)*100;
	nt=linspace(1,100,length(nFpp))';
	[p,ps] = polyfit(nt, nFpp, 1);
	R2 = (1 - ps.normr^2 / norm(nFpp-mean(nFpp))^2);
	
	
	if strcmp(plotflag,'plot') && length(h)>3
		axes(h(4));
		plot(nt,nFpp,'Color',c2)
		plot(nt, polyval([p(1) p(2) ],nt),'--','Color',c)
		legend off
	end
	
	
	%S
	ix = find(strcmp(parNames,'Sn'));
	if ix
		S = p(1);
		parValues(ix) = S;
		parValues_conf(ix) = R2;
		if strcmp(plotflag,'plot') && length(h)>3
			axes(h(4));
			text(40,90-plot_ix*5,sprintf('S=%.2f (%.2f)',S,R2),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		end
	end
	%Sintercept
	ix = find(strcmp(parNames,'Sintercept'));
	if ix
		% equation is flow = p1*time + p2
		% to find the time intercept of the line at flow = 0
		% => time = -p2/p1
		Sintercept = (p(2)-1) / p(1);
		parValues(ix) = Sintercept;
		parValues_conf(ix) = R2;
		if strcmp(plotflag,'plot') && length(h)>3
			axes(h(4));
			plot(Sintercept,0,'o','Color',c)
			text(Sintercept,-plot_ix*5,sprintf('Sintercept=%.2f',Sintercept),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		end
	end
	% PARAMETER S2 (2nd order fit)
	ix = find(strcmp(parNames,'S2'));
	if ix
		S2=nan;
		if length(nt)>3 && length(nFpp)>3
			[p,ps] = polyfit(nt, nFpp, 1);
			R2 = (1 - ps.normr^2 / norm(nFpp-mean(nFpp))^2);
			
		end
		parValues(ix) = p(1);
		parValues_conf(ix) = R2;
		if strcmp(plotflag,'plot') && length(h)>3
			plot( nt,polyval(p,nt),':','Color',c)
			axes(h(4));
			text(80,90-plot_ix*5,sprintf('S2=%.2f (%.2f)',S2,R2),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		end
	end
	
	% PARAMETER Sexp (2nd order fit)
	ix = find(strcmp(parNames,'SE'));
	if ix
		SE=nan;
		if length(nt)>3 && length(nFpp)>3
			[curveE,gofE] = fit(nt, nFpp, 'exp1');
			SE = curveE.b;
		end
		parValues(ix) = SE;
		parValues_conf(ix) = gofE.rsquare;
		if strcmp(plotflag,'plot') && length(h)>3
			plot( nt,curveE.a.*exp(nt.*curveE.b),':','Color',c)
			axes(h(4));
			text(80,90-plot_ix*5,sprintf('SE=%.5f (%.2f)',SE,gofE.rsquare),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		end
	end
	% PARAMETER Sexp (2nd order fit)
	ix = find(strcmp(parNames,'SP'));
	if ix
		SP=nan;
		SP_R2=nan;
		if length(nt)>3 && length(nFpp)>3
			try
				[curveP,gofP] = fit(nt, fliplr(nFpp')', 'power1');
				SP = curveP.b;
				SP_R2=gofP.rsquare;
			catch e
				SP=nan;
				SP_R2=nan;
			end
		end
		parValues(ix) = SP;
		parValues_conf(ix) = SP_R2;
		if strcmp(plotflag,'plot') && length(h)>3
			plot( nt,  fliplr((curveP.a.*nt.^curveP.b)' )  ,'-.','Color',c)
			axes(h(4));
			text(80,100-plot_ix*5,sprintf('SP=%.2f (%.2f)',SP,SP_R2),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
			
			
		end
	end
	
	
	% PARAMETER Sexp (2nd order fit)
	ix = find(strcmp(parNames,'SPn'));
	if ix
		SP2=nan;
		SP2_R2=nan;
		if length(nt)>3 && length(Fpp)>3
			try
				nt2 = linspace(0,(length(Fpp)-1)/fs,length(Fpp))';
				[curveP2,gofP2] = fit(nt2(2:end), fliplr(Fpp(2:end)')', 'power1');
				SP2 = curveP2.b;
				SP2_R2=gofP2.rsquare;
			catch e
				SP2 = nan;
				SP2_R2=nan;
			end
		end
		parValues(ix) = SP2;
		parValues_conf(ix) = SP2_R2;
		if strcmp(plotflag,'plot') && length(h)>3
			axes(h(6));
			hold on
			plot(nt2,Fpp,'Color',c2)
			plot( nt2,  fliplr((curveP2.a.*nt2.^curveP2.b)' )  ,'-.','Color',c)
			title(sprintf('SP=%.2f (%.2f)',SP2,SP2_R2),'Color',c,'VerticalAlignment','top','HorizontalAlignment','center','FontSize',fsize)
		end
	end
	
	% PARAMETER Sexp (2nd order fit)
	ix = find(strcmp(parNames,'SEn'));
	if ix
		SE2=nan;
		SE2_R2=nan;
		if length(nt)>3 && length(Fpp)>3
			try
				nt2 = linspace(0,(length(Fpp)-1)/fs,length(Fpp))';
				[curveP2,gofP2] = fit(nt2(2:end), Fpp(2:end), 'exp1');
				SE2 = curveP2.b;
				SE2_R2=gofP2.rsquare;
			catch e
				SE2 = nan;
				SE2_R2=nan;
			end
		end
		parValues(ix) = SE2;
		parValues_conf(ix) = SE2_R2;
		if strcmp(plotflag,'plot') && length(h)>3
			axes(h(6));
			hold on
			plot(nt2,Fpp,'Color',c2)
			plot( nt2,  curveP2.a.*exp(nt2.*curveP2.b)  ,'--','Color',c)
			text(0,0.3,sprintf('SE2=%.2f (%.2f)',SE2,SE2_R2),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end
	
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TRS EV 5                                                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETER FITTED STRIGHT LINE
if ~isempty(find(strcmp(parNames,'Trs'))) % || ~isempty(find(strcmp(parNames,'ar'))) || ~isempty(find(strcmp(parNames,'a0')))
	%prepare signal
	vrEx = vr(length(Vin)+1:end);
	frEx = fr(length(Vin)+1:end);
		

	[Cs,Rs,gama,SPvf,SPvf_R2,v_spvf,f_spvf] = FV_trs_ev_5gd(vrEx,frEx,'');
	
	
	TrsA_5 = -1/Cs(1,1);	
	TrsB_5 = -1/Cs(1,2);
	f_A0_5 = polyval(Cs(:,1),1);
	f_B0_5 = polyval(Cs(:,2),1);
	EV5 = f_A0_5*TrsA_5;
	Rab5 = sqrt(mean(Rs.^2));
	gama_5 = gama;
	
		
	if strcmp(plotflag,'plot')
		axes(h(5));
		hold on
		plot(vrEx,frEx,[c2 '-'])
		X = [vrEx'  1+EV5];
		plot(X,polyval(Cs(:,1),X),[c '-']);
		plot(X,polyval(Cs(:,2),X),[c '--']);
		
		xlim([0 EV5+1.01])
		ylim([0 max(frEx)+0.02])
			
		grid on
		grid minor
	end
	
	
	
	%Trs;
	ix = find(strcmp(parNames,'Trs'));
	if ix
		parValues(ix) = TrsA_5;
		parValues_conf(ix) = Rs(1);
		if strcmp(plotflag,'plot')
			axes(h(5));
			[~,a0] = max(frEx);
			text(vrEx(a0),frEx(a0)+plot_ix*plot_offset_flow,sprintf('TrsA=%.2f TrsB_5=%.2f',TrsA_5,TrsB_5),'Color',c,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',fsize)
			%title(sprintf('d=%.2f A=%.2f R2=%.2f',d,A,R2))
		end
	end
		%Trs;
	ix = find(strcmp(parNames,'TrsB'));
	if ix
		parValues(ix) = TrsB_5;
		parValues_conf(ix) = Rs(2);		
	end
	%EV
	ix = find(strcmp(parNames,'EV'));
	if ix		
		parValues(ix) = EV5;
		parValues_conf(ix) = Rab5;
		if strcmp(plotflag,'plot')
			axes(h(5));
			plot(EV5+1,0,'o--','Color',c)
			text(EV5+1,0-plot_ix*plot_offset_flow,sprintf('EV=%.4f',EV5),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end
	%EV
	ix = find(strcmp(parNames,'fA0'));
	if ix		
		parValues(ix) = f_A0_5;
		parValues_conf(ix) = Rab5;
		if strcmp(plotflag,'plot')
			axes(h(5));
			plot(1,f_A0_5,'o--','Color',c)
			text(1,f_A0_5-plot_ix*plot_offset_flow,sprintf('f_{A0}=%.4f',f_A0_5),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end
	%fA0_5
	ix = find(strcmp(parNames,'fB0'));
	if ix		
		parValues(ix) = f_B0_5;
		parValues_conf(ix) = Rab5;
		if strcmp(plotflag,'plot')
			axes(h(5));
			plot(1,f_B0_5,'o--','Color',c)
			text(1,f_B0_5-plot_ix*plot_offset_flow,sprintf('f_{B0}=%.4f',f_B0_5),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end
	
		
	%SPvf
	ix = find(strcmp(parNames,'SPvf'));
	if ix
		parValues(ix) = SPvf;
		parValues_conf(ix) = SPvf_R2;
		if strcmp(plotflag,'plot')
			axes(h(5));			
			plot(v_spvf,f_spvf,'-.','Color',c)			
			text(0,0,sprintf('SPvf=%.4f',SPvf),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end
	
	%SPvf
	ix = find(strcmp(parNames,'Gama'));
	if ix
		parValues(ix) = gama;
		parValues_conf(ix) = Rab5;
		if strcmp(plotflag,'plot')
			axes(h(5));
			text(0.5,0,sprintf('gama=%.4f',gama_5),'Color',c,'VerticalAlignment','top','HorizontalAlignment','left','FontSize',fsize)
		end
	end

	
	
end





end




