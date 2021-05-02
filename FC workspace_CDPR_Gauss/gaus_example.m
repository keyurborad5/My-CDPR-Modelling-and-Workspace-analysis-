clc
W_aux=[-0.5946    0.5946   -0.5946    0.5946    0.5946   -0.5946    0.5946   -0.5946;
    0.4460   -0.4460   -0.4460    0.4460   -0.4460    0.4460    0.4460   -0.4460;
   -0.6690    0.6690   -0.6690    0.6690   -0.6690    0.6690   -0.6690    0.6690;
   -0.6690    0.4460         0    0.2230         0    0.2230   -0.6690    0.4460;
   0   -0.0743   0   -0.0743    0.6690   -0.5946    0.6690   -0.5946;
    0.5946   -0.4460    0   -0.1487   -0.4460    0.5946   -0.1487         0]
W_aux38=W_aux(:,[1,2,4,8,5,3,6,7])
f=0;    
all_rank_aux=rank(W_aux38(:,:));
    if rank(W_aux38(:,:))==6    % we egnore matrix which are not full rank that is 6
        count = 0;
        for i=1:n-1
            
            if abs(W_aux38(i,i))<=1e-6  % this case runs when your pivot point turns to be zero 
                
                for q= i+1:n
                    if abs(W_aux38(q,i))>=1e-6      %we swap this row with the rows below it if they are not zero
                        count = count+1;
                        temp = W_aux38(i,:);
                        W_aux38(i,:)=W_aux38(q,:);
                        
                        W_aux38(q,:)=temp;
                        break
                    end
%               W(:,:,j)=[W(:,:,j);W(i,:,j)];
%               W(i,:,j)=[];
                end
            end
            if count~=0| abs(W_aux38(i,i))>=1e-6
                   W_aux38(i,:)=W_aux38(i,:)/W_aux38(i,i);  %actual code for Gauss elimination
                   m=W_aux38(i+1:n,i)/W_aux38(i,i);
                   W_aux38(i+1:n,:)=W_aux38(i+1:n,:)-m*W_aux38(i,:)
            end
        end
        f=f+1;
        oneD(f,:) = W_aux38(n,6:8); %this stores the last rows last 3 columns terms which we need to eval
    end