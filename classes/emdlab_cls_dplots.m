classdef emdlab_cls_dplots < handle

    properties

        r1 = 1.5;
        r2 = 4;
        r3 = 6;
        g = 0.3;
        lam_color = [0.7,0.7,0.7];
        colors = struct('blue','b','red','r','green',[30,157,68]/255,'magneta',[128,0,128]/255);

        arrowHeadWidth = 0.4;
        arrowHeadHeight = 0.3;

        frameIndex = 0;
        imageResolution = 400;

        rotorPosition = 0;

    end

    properties (Dependent=true)
        rc
    end

    methods

        function setFrameResolution(obj, value)
            obj.imageResolution = value;
        end

        function setArrowHeadSize(obj, width, height)
            obj.arrowHeadWidth = width;
            obj.arrowHeadHeight = height;
        end

        function y = get.rc(obj)
            y = (obj.r3 - obj.r2 - obj.g/2)*0.4/2;
        end

        function addRoundStator(obj, isMoving)

            if nargin<2, isMoving=false; end                

            hold on;
            t = linspace(0,2*pi,200);
            p = fill(obj.r3*cos(t),obj.r3*sin(t),obj.lam_color);
            p.UserData.isMoving = isMoving;
            p = fill((obj.r2+obj.g/2)*cos(t),(obj.r2+obj.g/2)*sin(t),'w');
            p.UserData.isMoving = isMoving;
            axis off equal;

        end

        function addRoundRotor(obj, isMoving)

            if nargin<2, isMoving=false; end                

            hold on;
            t = linspace(0,2*pi,200);
            p = fill((obj.r2-obj.g/2)*cos(t),(obj.r2-obj.g/2)*sin(t),obj.lam_color);
            p.UserData.isMoving = isMoving;
            p = fill(obj.r1*cos(t),obj.r1*sin(t),'w');
            p.UserData.isMoving = isMoving;
            axis off equal;

        end

        function addStatorCoilArmIn(obj, ang, tag, clr)
            obj.addCoilArm(ang, tag, clr, 1, 1, 0);
        end

        function addStatorCoilArmOut(obj, ang, tag, clr)
            obj.addCoilArm(ang, tag, clr, 0, 1, 0);
        end

        function addRotorCoilArmIn(obj, ang, tag, clr)
            obj.addCoilArm(ang, tag, clr, 1, 0, 1);
        end

        function addRotorCoilArmOut(obj, ang, tag, clr)
            obj.addCoilArm(ang, tag, clr, 0, 0, 1);
        end

        function addCoilArm(obj, ang, tag, clr, inOutFlag, rsFlag, isMoving)

            hold on;
            ang = deg2rad(ang);
            t = linspace(0,2*pi,50);

            if rsFlag
                rtmp = obj.r2 + obj.g/2 + 1.5*obj.rc;
            else
                rtmp = obj.r2 - obj.g/2 - 1.5*obj.rc;
            end

            x0 = rtmp*cos(ang);
            y0 = rtmp*sin(ang);
            x = x0 + obj.rc*cos(t);
            y = y0 + obj.rc*sin(t);
            p = fill(x,y,'w', 'LineWidth',1, 'EdgeColor',clr);
            p.UserData.isMoving = isMoving;

            if inOutFlag
                x = x0 + obj.rc*[-1,1]/sqrt(2);
                y = y0 + obj.rc*[-1,1]/sqrt(2);
                p = plot(x,y,'color',clr, 'LineWidth',1.5);
                p.UserData.isMoving = isMoving;
                y = y0 + obj.rc*[1,-1]/sqrt(2);
                p = plot(x,y,'color',clr, 'LineWidth',1.5);
                p.UserData.isMoving = isMoving;
            else
                p = plot(x0,y0,'color',clr, 'Marker','o','MarkerFaceColor',clr,'MarkerSize',5);
                p.UserData.isMoving = isMoving;
            end

            u = [x0,y0];
            u = u/norm(u);
            if rsFlag
                x = x0 + 2.3*obj.rc*u(1);
                y = y0 + 2.3*obj.rc*u(2);
            else
                x = x0 - 2.3*obj.rc*u(1);
                y = y0 - 2.3*obj.rc*u(2);
            end

            p = text(x,y,tag,'fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
            p.UserData.isMoving = isMoving;

        end

        function addAxis(obj, ang, tag, clr, scale, isMoving)

            if nargin<5, scale=1; end
            if nargin<6, isMoving=0; end
            hold on;
            ang = deg2rad(ang);

            x1 = scale * obj.r2 * cos(ang);
            y1 = scale * obj.r2 * sin(ang);

            u = [x1,y1];
            u = u/norm(u);
            n = [-u(2),u(1)];

            x2 = x1 - obj.arrowHeadWidth*0.5*u(1) + obj.arrowHeadHeight*0.5*n(1);
            y2 = y1 - obj.arrowHeadWidth*0.5*u(2) + obj.arrowHeadHeight*0.5*n(2);

            x3 = x1 - obj.arrowHeadWidth*0.5*u(1) - obj.arrowHeadHeight*0.5*n(1);
            y3 = y1 - obj.arrowHeadWidth*0.5*u(2) - obj.arrowHeadHeight*0.5*n(2);

            p = patch('Faces',[1,2,3,4,2], 'Vertices',[0,0;x1,y1;x2,y2;x3,y3],'linewidth',1,'edgecolor', clr ,'facecolor',clr);
            p.UserData.isMoving = isMoving;

            p = text(x1/2,y1/2,tag,'fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle');
            p.UserData.isMoving = isMoving;

        end

        function exportFrame(obj)

            obj.frameIndex = obj.frameIndex + 1;
            exportgraphics(gcf,'frame-'+string(obj.frameIndex)+'.png','Resolution',obj.imageResolution);

        end

        function rotate(obj, rotAngle)

            rotAngle = deg2rad(rotAngle);
            obj.rotorPosition = obj.rotorPosition + rotAngle;
            ax = gca;
            for i = 1:numel(ax.Children)
                p = ax.Children(i);
                if p.UserData.isMoving
                    if isa(p,'matlab.graphics.primitive.Text') 
                    x_new = p.Position(1)*cos(rotAngle) - p.Position(2)*sin(rotAngle);
                    y_new = p.Position(1)*sin(rotAngle) + p.Position(2)*cos(rotAngle);
                    p.Position(1) = x_new;
                    p.Position(2) = y_new;
                    elseif isa(p,'matlab.graphics.chart.primitive.Quiver')
                        x_new = p.XData*cos(rotAngle) - p.YData*sin(rotAngle);
                    y_new = p.XData*sin(rotAngle) + p.YData*cos(rotAngle);
                    p.XData = x_new;
                    p.YData = y_new;

                    x_new = p.UData*cos(rotAngle) - p.VData*sin(rotAngle);
                    y_new = p.UData*sin(rotAngle) + p.VData*cos(rotAngle);
                    p.UData = x_new;
                    p.VData = y_new;

                    else
                        x_new = p.XData*cos(rotAngle) - p.YData*sin(rotAngle);
                    y_new = p.XData*sin(rotAngle) + p.YData*cos(rotAngle);
                    p.XData = x_new;
                    p.YData = y_new;
                    end
                end
            end

        end

        function addAngle(obj, t1, t2, r, tag)

            hold on;
            t1 = deg2rad(t1);
            t2 = deg2rad(t2);
            t = linspace(t1, t2, 21);

            plot(r*cos(t),r*sin(t),'k', 'LineWidth',1.5);
            text(1.2*r*cos(t(11)),1.2*r*sin(t(11)),tag,'fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle');

        end

        function addRedAngle(obj, t1, t2, r, tag)

            hold on;
            t1 = deg2rad(t1);
            t2 = deg2rad(t2);
            t = linspace(t1, t2, 21);

            plot(r*cos(t),r*sin(t),'r', 'LineWidth',3);
            text(1.1*r*cos(t(11)),1.1*r*sin(t(11)),tag,'fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle', 'Color','r');

        end

        function addDashedCircle(obj, r)

            hold on;
            t = linspace(0,2*pi,100);
            p = plot(r*cos(t), r*sin(t), '--', 'Color', 'k');
            p.UserData.isMoving = false;
            axis off equal;

        end


        function addDashedLine(obj, r1, r2, ang)

            hold on;
            ang = deg2rad(ang);
            r = linspace(r1,r2,2);
            p = plot(r*cos(ang), r*sin(ang), '--', 'Color', 'k');
            p.UserData.isMoving = false;
            axis off equal;

        end

        function addCoil(obj, r, ang, tag, clr)

            if nargin<5, clr='k'; end

            hold on;
            ang = deg2rad(ang);
            width = 1;
            t = linspace(0,3.5*2*pi,500);
            x = cos(t)*width/4-t*width/14;
            y = sin(t)*width/2;
            x = [x(1),x,x(end)];
            y = [-width,y,-width];

            x = x + (max(x)-min(x))/2 + r - width/4;
            [x,y] = emdlab_g2d_rotatePointsXY(x,y,ang);
            plot(x,y, 'LineWidth',1.5, 'Color',clr)

            t = linspace(0,2*pi,10);
            fill(x(1) + 0.15*cos(t), y(1) + 0.15*sin(t), 'w', 'LineWidth',1.5, 'EdgeColor',clr)
            fill(x(end) + 0.15*cos(t), y(end) + 0.15*sin(t), 'w', 'LineWidth',1.5, 'EdgeColor',clr)

            u = [r*cos(ang), r*sin(ang)];
            n = [-u(2),u(1)]; n = n/norm(n);
            u = u + 1*width*n;
            text(u(1),u(2),tag,'fontsize',14,'HorizontalAlignment','center','VerticalAlignment','middle','Color',clr);


        end

        function addBgVectors(obj,scale,ang,clr,isMoving)

            hold on;
            ang = deg2rad(ang);
            t = linspace(0,2*pi,100);
            x = obj.r2*cos(t);
            y = obj.r2*sin(t);
            Bx = scale*cos(t-ang).*cos(t);
            By = scale*cos(t-ang).*sin(t);
            p = quiver(x,y,Bx,By,scale,'Color',clr);
            p.UserData.isMoving = isMoving;

        end


        function addBgVectors1(obj,scale,ang,clr,isMoving)

            hold on;
            ang = deg2rad(ang);
            t = linspace(0,2*pi,100);
            x = obj.r2*cos(t);
            y = obj.r2*sin(t);
            Bx = scale*cos(t-ang).*cos(t);
            By = scale*cos(t-ang).*sin(t);
            p = quiver(x,y,Bx,By,scale,'Color',clr);
            p.UserData.isMoving = isMoving;
            p = plot(x+Bx,y+By,'Color',clr);
            p.UserData.isMoving = isMoving;

        end

    end

end