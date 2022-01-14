clear
clc
% Change this line to your directory
Datas1 = load('/Users/apple/Documents/MATLAB/dataset/M-G_3D_x2.csv');
T1 = Datas1(:,1).';
X1 = Datas1(:,2).';
Y1 = Datas1(:,3).';

PointsCluster = zeros(size(X1));
Nnew = 10;
Clusters.X(1,:) = X1(1:Nnew);
Clusters.Y(1,:) = Y1(1:Nnew);
Clusters.Active(1) = true;
Clusters.Time(1,:) = T1(1:Nnew);

PointsCluster(1:Nnew) = 1;

Th_for_Merge_Clusters = 2;
Th_for_Add_Point = 6;
Std0 = 0.1;
TimeForget = 3;
for K = Nnew+1:numel(X1)
    tic;
    %%  Search for clusters to be merged
    NumCluster = numel(Clusters.Active);
    NumCluster1 = NumCluster;
    for K1 = 1:NumCluster-1
        for K2 = K1+1:NumCluster
            if Clusters.Active(K1)&&Clusters.Active(K2)% both are active
                % find mean and std of clusters
                Xm1 = mean(Clusters.X(K1,:)); Xm2 = mean(Clusters.X(K2,:)); Ym1 = mean(Clusters.Y(K1,:)); Ym2 = mean(Clusters.Y(K2,:));
                if numel(Clusters.X(K1,:))>1
                    Xs1 = std(Clusters.X(K1,:)); Ys1 = std(Clusters.Y(K1,:));
                else
                    Xs1 = Std0; Ys1 = Std0;
                end
                if numel(Clusters.X(K2,:))>1
                    Xs2 = std(Clusters.X(K2,:)); Ys2 = std(Clusters.Y(K2,:));
                else
                    Xs2 = Std0; Ys2 = Std0;
                end
                d = (Xm1-Xm2).^2/(Xs1*Xs2)+(Ym1-Ym2).^2/(Ys1*Ys2);
                if d < Th_for_Merge_Clusters
                    Clusters.Active(K1) = false; Clusters.Active(K2) = false; % both clusters are deacivated
                    NumCluster1 = NumCluster1 + 1;
                    Clusters.Active(NumCluster1) = true;
                    Xtemp = [Clusters.X(K1,:) Clusters.X(K2,:)]; Ytemp = [Clusters.Y(K1,:) Clusters.Y(K2,:)]; Ttemp = [Clusters.Time(K1,:) Clusters.Time(K2,:)];
                    [Ttemp,IS] = sort(Ttemp); Xtemp = Xtemp(IS); Ytemp = Ytemp(IS);
                    if numel(Ttemp) < Nnew; IxFirst = 1; else; IxFirst = numel(Ttemp)-Nnew+1; end
                    Clusters.X(NumCluster1,:) = Xtemp(IxFirst:end);
                    Clusters.Y(NumCluster1,:) = Ytemp(IxFirst:end);
                    Clusters.Time(NumCluster1,:) = Ttemp(IxFirst:end);
                end
            end
        end
    end
    %% Add points to the clusters
    ActiveClusters = find(Clusters.Active);
    Mx = []; My = []; Sx = []; Sy = []; Tend = [];
    for K3 = 1:numel(ActiveClusters)
        XX = Clusters.X(ActiveClusters(K3),:); XX(isnan(XX)) = [];
        YY = Clusters.Y(ActiveClusters(K3),:); YY(isnan(YY)) = [];
        Mx(K3) = mean(XX);My(K3) = mean(YY); %#ok<SAGROW>
        if numel(XX) > 1; Sx(K3) = std(XX);Sy(K3) = std(YY); else; Sx(K3) = Std0; Sy(K3) = Std0; end %#ok<SAGROW>
        Ttemp = Clusters.Time(ActiveClusters(K3),:);
        Ttemp = sort(Ttemp(~isnan(Ttemp)));
        Tend(K3) = Ttemp(end); %#ok<SAGROW>
    end
    Sx(Sx>Std0) = Std0; Sy(Sy>Std0)=Std0;
    d = (((X1(K)-Mx)./Sx).^2+((Y1(K)-My)./Sy).^2).*abs((T1(K)-Tend)./TimeForget).^0.25;
    [dmin,Imin] = min(d);
    if dmin < Th_for_Add_Point % Add point to the cluster
        PointsCluster(K) = ActiveClusters(Imin);
        SSX(K) = Sx(Imin);
        SSY(K) = Sy(Imin);
        MMX(K) = Mx(Imin);
        MMY(K) = My(Imin);
        % update cluster
        Xtemp = [Clusters.X(ActiveClusters(Imin),:) X1(K)]; Ytemp = [Clusters.Y(ActiveClusters(Imin),:) Y1(K)]; Ttemp = [Clusters.Time(ActiveClusters(Imin),:) T1(K)];
        Xtemp(isnan(Xtemp)) = []; Ytemp(isnan(Ytemp)) = []; Ttemp(isnan(Ttemp)) = []; 
        [Ttemp,IS] = sort(Ttemp); Xtemp = Xtemp(IS); Ytemp = Ytemp(IS);
        if numel(Ttemp) > Nnew
            Ttemp = Ttemp(numel(Ytemp)-Nnew+1:end);
            Xtemp = Xtemp(numel(Ytemp)-Nnew+1:end);
            Ytemp = Ytemp(numel(Ytemp)-Nnew+1:end);
        elseif numel(Ttemp) < Nnew
            Ttemp = [Ttemp NaN(1,Nnew-numel(Ytemp))]; %#ok<AGROW>
            Xtemp = [Xtemp NaN(1,Nnew-numel(Ytemp))]; %#ok<AGROW>
            Ytemp = [Ytemp NaN(1,Nnew-numel(Ytemp))]; %#ok<AGROW>
        end    
        Clusters.X(ActiveClusters(Imin),:) = Xtemp;
        Clusters.Y(ActiveClusters(Imin),:) = Ytemp;
        Clusters.Time(ActiveClusters(Imin),:) = Ttemp;
    else % produce new cluster
        NumClusters = numel(Clusters.Active)+1; 
        Clusters.Active(NumClusters) = true;
        Clusters.X(NumClusters,:) = [X1(K) NaN(1,Nnew-1)];
        Clusters.Y(NumClusters,:) = [Y1(K) NaN(1,Nnew-1)];
        Clusters.Time(NumClusters,:) = [T1(K) NaN(1,Nnew-1)];
        PointsCluster(K) = NumClusters;
    end
    te(K)=toc;
end

NC = unique(PointsCluster);
for K = 1:numel(NC)
    IX = find(PointsCluster == NC(K));
    plot3(T1(IX),X1(IX),Y1(IX),'.','MarkerSize',0.25);
    %plot3(IX,X1(IX),Y1(IX),'.','MarkerSize',0.25);
    hold all;%pause
end
grid on
