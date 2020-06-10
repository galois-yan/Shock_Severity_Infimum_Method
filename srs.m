classdef srs
    %SRS Summary of this class goes here
    %   Constructor function: 'SRS=srs(acc, StaF, Q)';
    %   Input 'acc' must be a object.
    
    properties
        fn
        StaF
        Q
        t
        resp
    end
    properties (Dependent)
        maximaxSRS
        tPeak
    end
    
    methods
        function value=get.maximaxSRS(SRS)
            value=max(abs(SRS.resp),[],2);
        end
        function value=get.tPeak(SRS)
            [M,I]=max(abs(SRS.resp),[],2);
            value=SRS.t(I);
        end
        function SRS=srs(fn, StaF, Q, t, resp)
            
            SRS.fn=fn';
            SRS.StaF=StaF;
            SRS.Q = Q;
            SRS.t=t;
            SRS.resp=resp;
        end
        function plot(SRS, varargin)
            figure;
            plot(SRS.fn, SRS.maximaxSRS);
            if nargin==2
                hold on;
                plot(SRS.fn, max(SRS.resp,[],2));
                plot(SRS.fn, -min(SRS.resp,[],2));
            end
            grid;
            set(gca,'MinorGridLineStyle',':','GridLineStyle',':','XScale','log','YScale','log');
            ylabel('Peak Accel (m/sec^2)');
            xlabel('Natural Frequency (Hz)');
            xlim([SRS.fn(1),SRS.fn(end)]);
        end
        
        function contourf(SRS)
            [X,Y]=meshgrid(SRS.t,SRS.fn);
            figure;
            resp=SRS.resp; %#ok<*PROP>

            %cmin=log10(SRS.StaF);
            cmin=log10(min(SRS.maximaxSRS));
            respLog=log10(abs(resp));
            respLog(respLog<cmin)=cmin;
            contourf(X,Y,respLog,'LineStyle','none');
            %surf(X,Y,respLog,'LineStyle','none');
            set(gca, 'YScale', 'log');
            colormap jet;
            c=colorbar;
            c.TickLabels=10.^c.Ticks;
            
            hold on;
            plot(SRS.tPeak,SRS.fn,'w','LineWidth',1.5);
            caxis([cmin,log10(max(SRS.maximaxSRS))]);
            xlabel('Time (s)');
            ylabel('Natural Frequency (Hz)');
        end
       
        function SRS=svd(SRS0,k)
            % 'k' is the order in SSM method.
            [U,S,V] = svd(abs(SRS0.resp'));
            x=diag(S);
            cost=(norm(x(2:end))/norm(x))^2;
            disp(['Cost function is ', num2str(cost)]);
            
            resp=U(:,k)*S(k,k)*V(:,k)';
            SRS=SRS0;
            SRS.resp=resp';
        if nargout==0
            L=abs(20*log10(SRS0.maximaxSRS./SRS.maximaxSRS));
            
            SRS.plot;
            hold on;
            plot(SRS.fn, SRS0.maximaxSRS);
            yyaxis right;
            plot(SRS.fn,L,'--','LineWidth',1);
            set(gca,{'ycolor'},{'k'});
            ylabel('Ratio (dB)');
            
            lgd =legend('SSI','SRS','Margin');
            lgd.Location='northwest';
            annotation('textbox',[0.68,0.18,0.186,0.058],'String',['L = ',num2str(mean(L)), ' dB'],'FitBoxToText','on');
            
            figure;
            for i=1:6
                subplot(5,2,i);
                plot(SRS0.t, -U(:,i)/max(abs(U(:,i))));
                xlabel({'Time, s'});
                ylabel({'Acceleration, m/s^2'});
                ylim([-1, 1]);
                legend(['\bf u_{',num2str(i),'}/|| \bf u_{',num2str(i),'}||_\infty']);
                set(gca,'GridLineStyle','--','XGrid','on','YGrid','on');
            end
            subplot(5,2,1);
            legend('\bf u_{ssi}');
        end
        end
        
        function [aM, aN]=mdof(SRS, MI)
            %MI is the modal information matrix with original sign, where MI(:,1) is natural
            %frequency, MI(:,1) is participation factor.
            x=MI(:,2).*MI(:,3);
            [X,Y]=meshgrid(SRS.fn, SRS.t);
            [Xq,Yq]=meshgrid(MI(:,1), SRS.t);
            Mq = interp2(X,Y,SRS.resp',Xq,Yq);
            yM=Mq*x;
            aM=acc([SRS.t, yM]);
            disp(['Absolute maximum response by matrix M is ', num2str(max(abs(yM)))]);
            
            Nq = interp2(X,Y,abs(SRS.resp'),Xq,Yq);
            yN=Nq*abs(x);
            aN=acc([SRS.t, yN]);
            disp(['Absolute maximum response by matrix N is ', num2str(max(yN))]);
            
            n = interp1(SRS.fn, SRS.maximaxSRS, MI(:,1));
            yn=n'*abs(x);
            disp(['Absolute maximum response by maximaxSRS is ', num2str(yn)]);
            
            SSM = SRS.svd(1);
            n1 = interp1(SSM.fn, SSM.maximaxSRS, MI(:,1));
            yn1=n1'*abs(x);
            disp(['Absolute maximum response by n1 is ', num2str(yn1)]);
        end
       
        
    end
    
    methods (Static)
        function[a1,a2,b1,b2,b3,rd_a1,rd_a2,rd_b1,rd_b2,rd_b3]=...
                srs_coefficients(f,damp,dt)
            ialgorithm=2;
            %
            tpi=2*pi;
            %
            num_fn=max(size(f));
            %
            a1=zeros(num_fn,1);
            a2=zeros(num_fn,1);
            b1=zeros(num_fn,1);
            b2=zeros(num_fn,1);
            b3=zeros(num_fn,1);
            %
            rd_a1=zeros(num_fn,1);
            rd_a2=zeros(num_fn,1);
            rd_b1=zeros(num_fn,1);
            rd_b2=zeros(num_fn,1);
            rd_b3=zeros(num_fn,1);
            %
            num_damp=length(damp);
            %
            for j=1:num_fn
                %
                omega=(tpi*f(j));
                %
                if(num_damp==1)
                    ddd=damp;
                else
                    ddd=damp(j);
                end
                %
                omegad=(omega*sqrt(1.-ddd^2));
                %
                cosd=(cos(omegad*dt));
                sind=(sin(omegad*dt));
                domegadt=(ddd*omega*dt);
                %
                rd_a1(j)=2.*exp(-domegadt)*cosd;
                rd_a2(j)=-exp(-2.*domegadt);
                rd_b1(j)=0.;
                rd_b2(j)=-(dt/omegad)*exp(-domegadt)*sind;
                rd_b3(j)=0;
                %
                if(ialgorithm==1)
                    %
                    a1(j)=(2.*exp(-domegadt)*cosd);
                    a2(j)=(-exp(-2.*domegadt));
                    b1(j)=(2.*domegadt);
                    b2(j)=(omega*dt*exp(-domegadt));
                    b2(j)=b2(j)*(( (omega/omegad)*(1-2*(ddd^2))*sind -2*ddd*cosd ));
                    b3(j)=0.;
                    %
                else
                    %
                    E=0;
                    K=0;
                    C=0;
                    S=0;
                    Sp=0;
                    %
                    E=(exp(-ddd*omega*dt));
                    K=(omegad*dt);
                    C=(E*cos(K));
                    S=(E*sin(K));
                    %
                    Sp=S/K;
                    %
                    a1(j)=(2*C);
                    a2(j)=(-(E^2));
                    %
                    b1(j)=(1.-Sp);
                    b2(j)=(2.*(Sp-C));
                    b3(j)=((E^2)-Sp);
                end
                %
            end
        end
        
    end
end
