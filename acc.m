classdef acc < timeseries
    %ACCELERATION Summary of this class goes here
    %   Constructor Method: acc(Sign); where 'Sign' contain some columns input.
    
    properties
        Sf  % Current sample rate.
    end
    
    methods
        function Acc=acc(Sign)
            Acc@timeseries(Sign(:,2:end), Sign(:,1)-Sign(1,1));
            Acc.Sf = (size(Sign,1)-1)/(Sign(end,1)-Sign(1,1));
        end
        %------------------------------------------------------------------
        function Acc=resample1(Acc, nSamples)
            t=linspace(Acc.Time(1),Acc.Time(end),nSamples)';
            y=Acc.Data;
            while 1
                try
                    y=resample(y, nSamples, length(y));
                    break;
                catch ME
                    y=y(1:2:end);
                end
            end
            Name=Acc.Name;
            Acc=acc([t,y]);
            Acc.Name=Name;
        end
        %------------------------------------------------------------------
        function SRS = srs(Acc, StaF, Q)
            t=Acc.Time;
            y=Acc.Data;
            tmx=max(t);
            tmi=min(t);
            nnn = length(t);
            dt=(tmx-tmi)/(nnn-1);
            sr=1./dt;
            
            fn(1)=StaF;
            if fn(1)>sr/30.
                fn(1)=sr/30.;
            end
            damp=1./(2.*Q);
            j=1;
            while(1)
                if (fn(j) > sr/8.)
                    break
                end
                fn(j+1)=fn(1)*(2. ^ (j*(1./24.)));
                j=j+1;
            end
            %
            [a1,a2,b1,b2,b3]=srs.srs_coefficients(fn,damp,dt);
            
            tmax=(tmx-tmi) + 1./fn(1);
            
            limit = round( tmax/dt );
            tt=[t',tmx+dt:dt:tmax];
            yy=[y',zeros(size(y,2),limit-nnn)];
            
            ns=max(size(fn));
            %
            resp=zeros(length(fn),length(t));
            for j=1:ns
                %
                forward=[ b1(j),  b2(j),  b3(j) ];
                back   =[     1, -a1(j), -a2(j) ];
                %
                resp(j,:)=filter(forward,back,y',[],2);
                %
            end
            
            SRS=srs(fn, StaF, Q, t, resp);
            if nargout==0
                SRS.plot;
            end
        end

        %------------------------------------------------------------------
        function plot(Acc)
            figure;
            plot(Acc.Time, Acc.Data);
            
            xlabel({'Time, s'});
            ylabel({'Acceleration, m/s^2'});
            set(gca,'GridLineStyle','--','XGrid','on','YGrid','on');
        end
        %------------------------------------------------------------------
        function [f, P1, P2, Y]=fft(Acc)
            Y=fft(Acc.Data);
            L=length(Acc.Data);
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1,:);
            P1(2:end-1) = 2*P1(2:end-1);
            f = Acc.Sf*(0:(L/2))/L;
            if nargout==0
                figure;
                plot(f,P1);
                title('Single-Sided Amplitude Spectrum');
                xlabel('\xi (Hz)');
                ylabel('$$|\hat{f}(\xi)|$$','Interpreter','Latex');
            end
        end
        %------------------------------------------------------------------
        
    end
    
end
