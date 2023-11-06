function out = KMC(X,Adj,rates,L,tmax,imax,plt,tmin)
%KMC Kinetic Monte Carlo (=Gillespie) algorithm for square lattice
%
%  IN:  X            x0 (N x nvars)
%       Adj          adjacency matrix (N x N)
%       rates.M1     conversion rates i->j (nvars x nvars)
%            .M2                     ik -> jk (nvars x nvars x nvars)
%       L            side of square lattice
%       tmax         max simulation time
%       imax         max # iterations
%       tmin         when to start recording output
%       plt.di       plot update interval
%          .clr          colours
%
%  OUT: out     matrix with [time, location, rind, tostate] for each reaction
%
% notes:
%   - taking averages of x,xy over time requires weighting 
%        each point i with Δti/ΣΔti (as dt/T in Integral[f(t)dt]/T)
%
% Bert Wuyts 26/06/2020

    if ~isempty(plt)
        figure('units','pixels','position',[50 50 500 500],'color','white');
        cmap = [0 0 0; cell2mat(plt.clr)];
    end

    %INITIALISATION
    N=size(X,1);
    M1=rates.M1;
    M2=rates.M2;
    nvars=size(M1,1);
    locs=(1:N)'; %coordinates
    nX=full(sum(X)); %initial type counts 
    
        %store neighbours (use cell allow generalisation to ntwks)
    nbrs=cell(N,1);
    for i=1:N
        nbrs{i}=locs(Adj(i,:));
    end
    
         %initialise output vars 
    out=nan(imax,4);  
    
        %count pairs at t=0
    IxORz=zeros(N,nvars); %surroundings of each cell with X or Z
    [i2,~,k2]=ind2sub(size(M2),find(M2));
    XZvals=unique([i2',k2']);
    for i=XZvals
        IxORz(:,i)=sum(Adj(:,X(:,i)),2); 
    end
    nXY=X'*IxORz; %surroundings totals (nvars x nvars)
    
    %initialise reaction rates
    R1=nX'.*M1;
    R2=permute(permute(M2,[1,3,2]).*nXY,[1,3,2]);
        
    %GILLESPIE ALGORITHM (= KINETIC MONTE CARLO)
    %(=one reaction at a time=>only update relevant subsets)
    i=1;t=0;ii=0;
    
    while (t<tmax) && (i<imax)     
        %determine next reaction
        R=[R1(:);R2(:)];
        cumR=cumsum(R);
        U=rand;
        rind=find(cumR>U*cumR(end),1,'first'); 
        if isempty(rind)
            break %if no reactions left, stop
        elseif rind==1
            dU=cumR(end)*U/R(rind); 
        else %dU is how much U exceeds thrh of reaction rind, scaled (to determine next cell)
            dU=(cumR(end)*U-cumR(rind-1))/R(rind);
        end
        
        %update iteration and time
        i=i+1; 
        dt=-log(rand)/cumR(end);%dt ~ Exp[cumR(end)]: use inv tr of unif rand
        t=t+dt;
        
        %find location 
        if rind < nvars^2+1 
            %cell for X->Y
            [rindi,rindj]=ind2sub(size(M1),rind);
            locsi=locs(X(:,rindi));
            loc=locsi(floor(dU*nX(rindi))+1);
        else 
            %cell for XZ->YZ
            [rindi,rindj,rindk]=ind2sub(size(M2),rind-nvars^2);
            loc=find(cumsum(X(:,rindi).*IxORz(:,rindk))>dU*nXY(rindi,rindk),1,'first'); %!bottleneck! 
        end
        
        %output
        if t>=tmin
            ii=ii+1;
            if ii==1
                out(1:N,1)=t;
                out(1:N,2)=locs;
                out(1:N,3)=NaN;
                out(1:N,4)=sum(X*diag(1:nvars),2);
                ii=N+1;
            end
            out(ii,1)=t;
            out(ii,2)=loc;
            out(ii,3)=rind;
            out(ii,4)=rindj;
        end
        
        %update surroundings
        binXZ=ismember([rindi,rindj],XZvals);
        if any(binXZ)
            %store subset of nXY that will change
            nXYsub0=X([loc;nbrs{loc}],:)'*IxORz([loc;nbrs{loc}],:);
            %update IxORz
            if binXZ(1)
                IxORz(nbrs{loc},rindi)=IxORz(nbrs{loc},rindi)-1;
            end

            if binXZ(2)
                IxORz(nbrs{loc},rindj)=IxORz(nbrs{loc},rindj)+1;
            end
        end
                
        %update selected cell, nX  
        X(loc,rindi)=false;X(loc,rindj)=true;
        nX(rindi)=nX(rindi)-1;
        nX(rindj)=nX(rindj)+1;
        
        %update R1
        R1(rindi,:)=nX(rindi).*M1(rindi,:);
        R1(rindj,:)=nX(rindj).*M1(rindj,:);
        
        %update nXY and R2
        if any(binXZ)
            nXYsub1=X([loc;nbrs{loc}],:)'*IxORz([loc;nbrs{loc}],:);
            dnXY=nXYsub1-nXYsub0;%double check if this and next step really necessary
            nXY=nXY+dnXY;
            [indd1,indd2]=ind2sub([nvars,nvars],find(dnXY));
            R2(indd1,:,indd2)=permute(permute(M2(indd1,:,indd2),[1,3,2]).*nXY(indd1,indd2),[1,3,2]);
        end
             
        %plot
        if ~isempty(plt) && mod(i,plt.di)==0 %update plot
            imagesc(reshape(sum(X*diag(1:nvars),2),L,L));
            colormap(gca,cmap);caxis([0,nvars]);
            axis square; axis off;
            box on;
            title(sprintf('t=%.6f/%u',t,tmax));
            drawnow;%waitforbuttonpress;%
        end
    end
    
    out=out(1:ii,:); %remove superfluous rows
end
