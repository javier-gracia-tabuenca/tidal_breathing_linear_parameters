function [Cs,Rs,gama,SPvf,SPvf_R2,v_spvf,f_spvf] = FV_trs_ev_5(v,f,varargin)
% Calculates ar and a0 as seen in [1] from a expiration ciclein a
% flow-volume loop. 
% double line fit ??

plotflag='';
n = 0;
ndownsamp = 0;
while n < length(varargin)
    n = n + 1;
    if strcmp(varargin{n}, 'plot')
        plotflag = 'plot';
% 	elseif strcmp(varargin{n}, 'weightA0')
% 		n = n + 1;
% 		w_a0 = varargin{n};
% 	elseif strcmp(varargin{n}, 'weightAr')
% 		n = n + 1;
% 		w_ar = varargin{n};
    end
end




%% resample in fv domain 
v0=v;
f0=f;

v = linspace(0,1,500);
f = interp1(v0,f0,v);
	
[~,ix_low] = max(f);
% ix_low = ix_low+100;
[~,ix_high] = min(abs(v-0.95));

vv = v(ix_low:ix_high);
ff = f(ix_low:ix_high);

if length(vv) < 3
	Cs=[0 0 ;0 0 ];
	Rs = [0 0 ];
	gama=0;
	SPvf=0;
	SPvf_R2=0;
	return
end

vv = vv-vv(1);
ff =  fliplr(ff);
ff = ff-ff(1);


vv = vv(2:end);
ff = ff(2:end);

[curveP2,gofP2] = fit(vv',ff', 'power1');
SPvf = curveP2.b;
SPvf_R2=gofP2.rsquare;

ffe = (curveP2.a.*vv.^curveP2.b); 

v_spvf = vv+v(ix_low);
f_spvf = fliplr(ffe)+f(ix_high);

%FIT " LINES
%pars
vv_limL = 0.15;%perone
vv_limH = 0.85;%perone
ff_limL = 0;%perone
ff_limH = 1;%perone

lvv=length(vv);
lff=length(ff);

X = linspace(vv(1),vv(end),lvv);
Y = linspace(ff(1),ff(end),lff);

xdata = [vv ; ff];
ydata = ff;

lb = [X(round(lvv*vv_limL)+1) Y(round(lff*ff_limL)+1)];
ub = [X(round(lvv*vv_limH)) Y(round(lff*ff_limH))];

x0 = [ X(round(lvv/2)) Y(round(lff/2)) ];

%Run lsq funation 
[x,resnorm,res]  = lsqcurvefit(@lsqfun,x0,xdata,ydata,lb,ub,optimoptions('lsqcurvefit','Display','off'));
% R squared
R2 = 1 - ( sum(res.^2) / sum((ff-mean(ff)).^2)); 

%repeat for solution 
[F,bl,br,R2l,R2r] = lsqfun(x,xdata);

al =  x(2) - x(1)*bl;
ar =  x(2) - x(1)*br;

c1 = [ bl al ];
c2 = [ br ar ];
r1 = R2l;
r2 = R2r;


Cs0 = [c1 ; c2]';
Rs = [r1 r2];
%gama =  angle between lines
gama = atan(Cs0(1,1))  - atan(Cs0(1,2));


%% Convert lines to the original fv loop
v00 = v(ix_high);
f00 = f(ix_high);
%
Cs = Cs0;
Cs(2,:) = Cs(2,:) + f00;
for i=1:2
	Cs(2,i) = polyval(Cs(:,i),v00) ;
end
Cs(1,:) = -Cs(1,:);



% 
% figure
% plot(r1)
% hold on 
% plot(r2)
% plot(m,'.-')
% [a,ix]=max(m);
% plot(ix,a,'o')
 
if strcmp(plotflag,'plot')
	
% 
% figure
% hold on
% plot(vv,ff)
% plot(X,polyval(Cs0(:,1),X),'b')
% plot(X,polyval(Cs0(:,2),X),'c')
% 
% 
% 
% title(sprintf('b= %f1.3 R= %f1.3',SPvf,SPvf_R2))
% 
% % 
% % figure
% % plot(v,f,'.')
% % hold on 
% % plot(vv,ff,'.')
% 
% grid on 
% grid minor
% 
figure
hold on

%contour
res = 3;

X2 = X(round(lvv*vv_limL)+1: res : round(lvv*vv_limH));
Y2 = Y(round(lff*ff_limL)+1: res : round(lff*ff_limH));


for ivv = 1:length(X2)
	for iff = 1:length(Y2)
		xx = [X2(ivv) Y2(iff)];
		F = lsqfun(xx,xdata);
		OP(ivv,iff) = sum( (F-ydata).^2);
	end
end

%norm
OP = OP./max(max(OP));
OP = flip(OP');
OP = flip(OP');
OP = flip(OP');



cm = flip(gray);

cm = cm - 0.2;
cm(find(cm<0))=0;

rectangle('Position',[v(ix_low)+lb(1) f(ix_high)+lb(2) ub(1)-lb(1) ub(2)-lb(2)],'FaceColor',[ 1 1 1]*0.97,'EdgeColor',[ 1 1 1]*0.99)
colormap(cm);
contour(X2+v(ix_low),Y2+f00,OP.^0.5,[0.035 0.1:0.1:1],'-')



%calculate intersection
x_intersect = fzero(@(x) polyval(Cs(:,1)-Cs(:,2),x),3);
y_intersect = polyval(Cs(:,1),x_intersect);


plot(v0,f0,'Color',[1 1 1]*0.8,'LineWidth',3)

%left
vl = [0 x_intersect];
plot(vl,polyval(Cs(:,2),vl),'k-');

vl = [x_intersect 1.4];
plot(vl,polyval(Cs(:,1),vl),'k-');


plot([-0.1 1.1],[f(ix_high) f(ix_high)],'k--')
plot([v(ix_low) v(ix_low)],[0 1.4],'k--')
plot([v(ix_high) v(ix_high)],[0 1.4],'k--')


xlim([-0.1 1.1])
ylim([0 f(ix_low)*1.2])


xticks([0 v(ix_low) v(ix_high) 1])
xticklabels({'0' 'V_{PTEF}' '0.97' ''})
yticks([])


fs=10;
drawbrace([v(ix_low)+lb(1) f(ix_high)], [v(ix_low) f(ix_high)], 10, 'Color', 'k');
text(v(ix_low)+lb(1)/2, f(ix_high)-0.07,'10%','HorizontalAlignment','center','FontSize',fs);

drawbrace([v(ix_low)+ub(1) f(ix_high)], [v(ix_low)+0.1 f(ix_high)], 10, 'Color', 'k');
text(v(ix_low)+(lb(1)+ub(1))/2, f(ix_high)-0.07,'80%','HorizontalAlignment','center','FontSize',fs);

drawbrace([v(ix_high) f(ix_high)], [v(ix_low)+ub(1) f(ix_high)], 10, 'Color', 'k');
text(v(ix_high)-(v(ix_high)-(v(ix_low)+ub(1)))/2, f(ix_high)-0.07,'10%','HorizontalAlignment','center','FontSize',fs);



xlabel('Vol [l]')
ylabel('flow [l/s]')
%%





end

end

function [F,bl,br,R2l,R2r] = lsqfun(x,xdata)
	
	vv = xdata(1,:);
	ff = xdata(2,:);

	v = x(1);
	f = x(2);
	
	Vl = vv(find(vv <= v)) - v;
	Vr = vv(find(vv  > v)) - v;
	%
	Fl = ff(find(vv <= v)) - f;
	Fr = ff(find(vv  > v)) - f;
		
	%fit left
	bl = Vl'\Fl';% x = A\B Solve systems of linear equations Ax = B for x
	fitFl = bl.*Vl+f;
		%fit rigth
	br = Vr'\Fr';% x = A\B Solve systems of linear equations Ax = B for x
	fitFr = br.*Vr+f;
	
	F = [fitFl fitFr];
	
	R2l = 1 - ( sum((Fl-(fitFl-f)).^2) / sum((Fl-mean(Fl)).^2));
	R2r = 1 - ( sum((Fr-(fitFr-f)).^2) / sum((Fr-mean(Fr)).^2));
	
% 	figure
% 	hold on 
% 	plot(vv,ff)
% 	X = linspace(vv(1),vv(end),length(vv));
% 	plot(X,F)
	
	
end

















